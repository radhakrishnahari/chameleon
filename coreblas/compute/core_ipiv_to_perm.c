/**
 *
 * @file core_ipiv_to_perm.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ipiv_to_perm CPU kernel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matteo Marcos
 * @date 2025-03-24
 */
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * The idea here is to generate a permutation from the sequence of
 * pivot.  To avoid storing one whole column at each step, we keep
 * track of two vectors of nb elements, the first one contains the
 * permutation of the first nb elements, and the second one contains
 * the inverse permutation of those same elements.
 *
 * Lets have i the element to pivot with ip. ipiv[i] = ip;
 * We set i_1 as such invp[ i_1  ] = i
 *  and  ip_1 as such invp[ ip_1 ] = ip
 *
 * At each step we want to:
 *   - swap perm[i] and perm[ip]
 *   - set invp[i_1] to ip
 *   - set invp[ip_1] to i
 *
 *******************************************************************************
 *
 * @param[in] m0
 *          The base index for all values in ipiv, perm and invp. m0 >= 0.
 *
 * @param[in] m
 *          The number of elements in perm and invp. m >= 0.
 *
 * @param[in] k
 *          The number of elements in ipiv. k >= 0.
 *
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done.
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done.
 *
 * @param[in] ipiv
 *          The pivot array of size n. This is a (m0+1)-based indices array to follow
 *          the Fortran standard.
 *
 * @param[out] perm
 *          The permutation array of the destination row indices (m0-based) of the [1,n] set of rows.
 *
 * @param[out] invp
 *          The permutation array of the origin row indices (m0-based) of the [1,n] set of rows.
 *
 */
void CORE_ipiv_to_perm( int m0, int m, int k, int K1, int K2, int *ipiv, int *perm, int *invp )
{
    int i, j, ip;
    int i_1, ip_1;

    for(i=0; i < m; i++) {
        perm[i] = i + m0;
        invp[i] = i + m0;
    }

    for(i = 0; i < k; i++) {
        if ( ( i + m0 < K1 ) || ( i + m0 > K2 ) ) {
            continue;
        }
        ip = ipiv[i]-1;
        assert( ip - m0 >= i );

        if ( ip - m0 > i ) {

            i_1 = perm[i];

            if (ip-m0 < m) {
                ip_1 = perm[ip-m0];
                perm[ip-m0] = i_1;
            } else {
                ip_1 = ip;
                for(j=0; j < m; j++) {
                    if( invp[j] == ip ) {
                        ip_1 = j + m0;
                        break;
                    }
                }
            }

            perm[i] = ip_1;
            i_1  -= m0;
            ip_1 -= m0;

            if (i_1  < m) invp[i_1 ] = ip;
            if (ip_1 < m) invp[ip_1] = i + m0;
        }
    }
}
