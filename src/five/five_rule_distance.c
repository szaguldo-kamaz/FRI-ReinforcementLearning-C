/*
** Fuzzy Inference by Interpolation in Vague Environment toolbox for ANSI C
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** ANSI C / AVX port and FixRes/NoInf/NoNan/etc. optimization by David Vincze <david.vincze@webcode.hu>
** Various contributions by Daniel Palko <palko.daniel@uni-miskolc.hu>
** mex port by Sandor Zsuga (zsuga@iit.uni-miskolc.hu)
** Original FIVE MATLAB code by Szilveszter Kovacs <szkovacs@iit.uni-miskolc.hu>
**
**
**FIVERuleDist: Calculate the distances of the observation from the rule antecedents
**
**                            [RD]=FIVERuleDist(U,VE,R,X)
**
**          Where:
**          U: is the universe (a vector of discrete values in increasing order),
**          VE: is the generated vague environment (the scaled distances on the universe U)
**             If VE(k,1)=0, it is the primitive integral of the scaling function according
**               to the elements of U
**             If VE(k,1)=-1, it is the integral of the scaling function between
**               the neighboring elements of U
**          R: is the Rulebase
**          X: is the observation
**          In case of U,VE,X, the rows are the dimensions
**          RD: is the scaled distance of the observation from the rule antecedents
**               RD<0 denotes that RD=inf. This case abs(RD) is a distance
**               counted based on the number of inf values in VE between the
**               observation from the rule antecedents.
**
**          Rulebase:
**
**            R1: a1 a2 a3 ... am -> b
**            R2: a1 a2 a3 ... am -> b
**            R3: a1 a2 a3 ... am -> b
**             ...
**            Rn: a1 a2 a3 ... am -> b
**
**            If ai=NaN then the rule antecedent ai is not given
*/

///
///\file five_rule_distance.c
///

#include "FIVE.h"

#include "../inl/fast_sqrt.inl"
#include "../inl/eucledian.inl"

#include "../inl/fast_abs.inl"
#include "../inl/arr.inl"
#include "../inl/min.inl"
#include "../inl/sort.inl"

#include <string.h>


/// Calculate the distances of the observation from the rule antecedents
///\param frb five struct
///\param x observation point
///\return void
int five_rule_distance(struct FIVERB *frb, fri_float *x) {

    // cartpole 124763x
    // acrobot   65868x

    fri_float *diststmp;

#if defined(BUILD_AVX2) && defined(FIVE_NONAN) && defined(FIVE_NOINF)

    unsigned int minjs[frb->numofunivs];

    for (unsigned int k = 0; k < frb->numofunivs; k++) {
        minjs[k] = get_vag_abs_min_i_fixres(frb->uk[k], frb->univlength, x[k], frb->udivs[k]);
//    }

//    for (unsigned int k = 0; k < frb->numofunivs; k++) {

        diststmp = frb->frd_dists + k*4; // 4*64bit = 256bit (ymm)

        asm volatile (
                "mov      %3, %%cx                 \n"  // cx = frb->avx2_rbsize | arrsize
                "mov      %2, %%r9                 \n"  // &frb->rseqant_veval[k][0]
                "mov      %1, %%r10                \n"  // &diststmp[0]
                "vbroadcastsd %0, %%ymm1           \n"  // ymm1 = < vek[k][minj] vek[k][minj] vek[k][minj] vek[k][minj] >

// no need for abs, will be squared anyway
        "frdl:   vmovdqa  (%%r9), %%ymm0           \n"  // ymm0 = < frb->rseqant_veval[k][i] i+1 i+2 i+3 >
                "vsubpd   %%ymm0, %%ymm1, %%ymm0   \n"  // ymm0 = ymm1 - ymm0
                "vmulpd   %%ymm0, %%ymm0, %%ymm0   \n"  // ymm0 = ymm0 * ymm0 // these will be summed later
                "vmovdqa  %%ymm0, (%%r10)          \n"  // diststmp[0] = ymm0

                "add     $32, %%r9         \n"
                "add    $256, %%r10        \n"  // 256 - 4 * FIVE_MAX_NUM_OF_UNIVERSES * sizeof(double)
//                "dec      %%cx             \n"  // sub $1 faster??
                "sub      $1, %%cx         \n"
                "jne      frdl\n"

        : : "m" (frb->vek[k][minjs[k]]), "m" (diststmp), "m" (frb->rseqant_veval[k]), "m" (frb->avx2_rbsize) : "r9", "r10", "rcx", "ymm0", "ymm1", "memory" );

    }

#else

    // This is taken from vagdist! FRIRL does not use INF or NAN
    for (unsigned int k = 0; k < frb->numofunivs; k++) {

        diststmp = frb->frd_dists + k*4;  // 4*64bit = 256bit (ymm)

        unsigned int minj = get_vag_abs_min_i_fixres(frb->uk[k], frb->univlength, x[k], frb->udivs[k]);

        for (unsigned int r = 0; r < frb->numofrules; r += 4) {  // 4*64bit = 256bit (ymm)
            // no need for abs(), because it will be squared anyway

            diststmp[0] = (frb->vek[k][frb->rseqant_uindex[k][r+0]] - frb->vek[k][minj]);
            diststmp[1] = (frb->vek[k][frb->rseqant_uindex[k][r+1]] - frb->vek[k][minj]);
            diststmp[2] = (frb->vek[k][frb->rseqant_uindex[k][r+2]] - frb->vek[k][minj]);
            diststmp[3] = (frb->vek[k][frb->rseqant_uindex[k][r+3]] - frb->vek[k][minj]);

            diststmp += 4 * FIVE_MAX_NUM_OF_UNIVERSES;
        }
    }

#endif // AVX2 or not

    /*
     *    dists =    r1,d1 r2,d1 r3,d1 r4,d1 r1,d2 r2,d2 r3,d2 r4,d2 ...
     *
     *        ... ri,dj, ri+1,dj, ri+2,dj, ri+3,dj, ri+1,dj+1, ri+2,dj+1, ri+3,dj+1, ri+4,dj+1 ...
     *
     *        ... rn-7,dm-1 rn-6,dm-1 rn-5,dm-1 rn-4,dm-1 rn-3,dm rn-2,dm rn-1;dm rn,dm
     *    r=(1..i..n), n=umofrules; d=(1..j..m), m=numofunivs
     */

#ifdef BUILD_AVX2

    /*
     *    arrsize = numofrules / 4
     *    i = j = k = 1
     *    sum = 0
     *
     *    while i < arrsize
     *        sum += dists[j,    ..., j+3] * dists[j,     ..., j+3]
     *        sum += dists[j+4,  ..., j+7] * dists[j+4,   ..., j+7]
     *        ...
     *        sum += dists[j+60, ..., j+63] * dists[j+60, ..., j+63]
     *
     *        ruledists[k, ..., k+3] = sqrt(sum)
     *    inc i
     *    j += 64
     *    k += 4
     */

    unsigned long int exacthit[4] __attribute__ ((aligned (32))) = { 1, 2, 3, 4 };
    unsigned int cnt=0;
    unsigned char rp;
    unsigned int rno;

    asm volatile (
        "mov %0,    %%rax;"  // &ruledists[0]
        "mov %1,    %%rbx;"  // &dists[0]
        "mov %2,    %%cx;"   // frb->avx2_rbsize | arrsize

#ifdef FRIRL_FAST
        "vpxor %%ymm0, %%ymm0, %%ymm0 \n"
#endif
        "l%=: ;"

        // FIVE_MAX_NUM_OF_UNIVERSES = 8 - unrolled loop for 4 rules at once - rule n...n+3
        "vmovapd (%%rbx), %%ymm3;"       // ymm3 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmovapd (%%rbx), %%ymm1;"       // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm3;" // ymm3 = ymm1 * ymm1

//        "vmovapd 32(%%rbx), %%ymm1;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
        "vmovapd 32(%%rbx), %%ymm2;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
        "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//        "vmovapd 64(%%rbx), %%ymm1;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
        "vmovapd 64(%%rbx), %%ymm2;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
        "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//        "vmovapd 96(%%rbx), %%ymm1;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
        "vmovapd 96(%%rbx), %%ymm2;"     // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
        "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//        "vmovapd 128(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
        "vmovapd 128(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
        "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//        "vmovapd 160(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
        "vmovapd 160(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
        "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//        "vmovapd 192(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
        "vmovapd 192(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
        "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2

//        "vmovapd 224(%%rbx), %%ymm1;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
        "vmovapd 224(%%rbx), %%ymm2;"    // ymm1 = <dists[i+3], dists[i+2], dist[i+1], dist[i]>
//        "vmulpd %%ymm1, %%ymm1, %%ymm2;" // ymm2 = ymm1 * ymm1
        "vaddpd %%ymm2, %%ymm3, %%ymm3;" // ymm3 += ymm2
        // unroll end

        "vsqrtpd %%ymm3,  %%ymm3;" // ymm3 = sqrt(ymm3)
        "vmovapd %%ymm3, (%%rax);" // frb->ruledists[i, ..., i+3] = ymm3
#ifdef FRIRL_FAST
    // exacthit - should be after sqrt() because other values in ymm3 can be useful values
        "vcmpeqpd %%ymm3, %%ymm0, %%ymm4  \n" // ymm3 van-e benne 0.0?
        "vptest %%ymm4, %%ymm4            \n" // csak flageket valtoztat (hogy az egesz ymm4 nulla-e)
        "jnz     frd_exact                \n" // ugrik ha valahol igaz, szoval ymm3 valahol 0.0 -> exact hit van
#endif
        "add $32,   %%rax;" // ruledists += 4
        "add $256,  %%rbx;" // dists += 32
        "sub $1,    %%cx;"  // i--, sub faster than dec??
        "jnz l%= ;"         // if arrsize != i, then goto l

        "jmp frd_ki \n"
#ifdef FRIRL_FAST
    "frd_exact: \n"
//    "vmovmskpd "
        "vmovapd %%ymm4, %3\n"
#endif
    "frd_ki: \n"
        "mov %%cx, %4 \n"

        : // no output, frb->ruledists directly updated
        : "m" (frb->ruledists), "m" (frb->frd_dists), "m" (frb->avx2_rbsize), "m" (exacthit), "m" (cnt)
        : "rax", "rbx", "rcx", "ymm0", "ymm2", "ymm3", "ymm4", "memory"
    );

#ifdef FRIRL_FAST
    // if loop exited before the last rule (probably an exact hit)
    // could be optimized a little, see and test vmovmskpd
    if (cnt > 0) {
        if (exacthit[0] > 0) {
            rp = 0;
        } else {
            if (exacthit[1] > 0) {
                rp = 1;
            } else {
                if (exacthit[2] > 0) {
                    rp = 2;
                } else {
                    if (exacthit[3] > 0) {
                        rp = 3;
                    }
                }
            }
        }

        rno = (frb->avx2_rbsize - cnt)*4 + rp;
        if (rno < frb->numofrules) {
            return rno;
        }
    }
#endif

    // debug info: a tobb 0.0 tavolsag ne tevesszen meg, azok az elozo iteraciokbol maradtak!
    // debug info: nem gond ha az ymm regben a dist[] input jon valami szemet es a vegeredmeny is szemet (nem 0.0), mert nincs figyelembe vev (rno < frb->numofrules)

#else // C code here

    diststmp = frb->frd_dists;
    for (unsigned int i = 0; i < frb->numofrules; i += 4) {
        for (unsigned int j = 0; j < 4; j++) {
            fri_float ret = 0.0;

            for (unsigned int k = 0; k < frb->numofunivs; k++) {
                fri_float dist = diststmp[k * 4 + j];
                ret += dist * dist;
            }

#ifdef FRIRL_FAST
            // if a direct hit is found then return with the rule number
            if ((ret == 0.0) && ((i+j) < frb->numofrules)) {
                return (i+j);
            }
#endif

            frb->ruledists[i + j] = fast_sqrt(ret);
        }
        diststmp += 4 * FIVE_MAX_NUM_OF_UNIVERSES;
    }

#endif

    return ~0;
}
