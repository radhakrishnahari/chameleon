/**
 *
 * @file zhmat2rpk.cc
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon tile to RAPACK format conversion
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-05-18
 * @precisions normal z -> c d s
 *
 */

#include "runtime_rpk.h"
#include "runtime_lrmem.h"
#include "chameleon/config.h"
#include "chameleon/constants.h"
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "coreblas.h"
#include <rapack.h>
#include <starpu.h>
#include "cham_tile_interface.h"


#if defined(CHAMELEON_USE_HMATOSS)
#include "runtime_hmat.h"
#include "coreblas/include/coreblas/hmat.h"
#include "coreblas/include/coreblas/hmat.hpp"
#include "coreblas/hmat-oss/src/data_types.hpp"
#include "coreblas/hmat-oss/src/tree.hpp"
#include "coreblas/hmat-oss/src/h_matrix.hpp"

using namespace hmat;

namespace {
    class CHAMELEON_ZHmat2Rpk_Proc : TreeProcedure<HMatrix<HMAT_CHAMELEON_Z_t>> {
    public:
        CHAM_tile_t *tileA;
        CHAM_tile_t *tileB;

        CHAMELEON_ZHmat2Rpk_Proc( CHAM_tile_t * tileA, CHAM_tile_t * tileB ) : tileA(tileA), tileB(tileB) {
            
        }

        virtual void visit( HMatrix<HMAT_CHAMELEON_Z_t>* node, const Visit order ) const {
            if ( !node->isLeaf() ) {
                fprintf(stderr, "error: hierarchical matrix conversion are not supported yet\n");
                exit(2);
            }
            assert( order == tree_leaf );
            if ( node->isFullMatrix() ) {
                FullMatrix<HMAT_CHAMELEON_Z_t> *fullmat = node->full();
                ScalarArray<HMAT_CHAMELEON_Z_t> *a_arr = &fullmat->data;
                const rpk_ctx_t *ctx = runtime_rpk_zctx_get();
                rpk_matrix_t *Bra = (rpk_matrix_t*)CHAM_tile_get_ptr( tileB );
                CHAMELEON_Complex64_t *u = (CHAMELEON_Complex64_t*)a_arr->ptr();
                assert( sizeof(HMAT_CHAMELEON_Z_t) == sizeof(CHAMELEON_Complex64_t) );

                /* Resize the rapack matrix to full rank */
                rpkx_zlrsze(
                    ctx,
                    0,
                    tileB->m, tileB->n,
                    Bra,
                    -1, -1,
                    -1
                );

                /* Set the content of the new rapack matrix */
                CORE_zlacpy(
                    ChamUpperLower,
                    tileA->m, tileA->n,
                    u, a_arr->lda,
                    (CHAMELEON_Complex64_t*)Bra->u, tileB->m
                );
            }
            else {
                assert( node->isRkMatrix() );
                int rank = node->rank();
                CHAMELEON_Complex64_t *a_u;
                CHAMELEON_Complex64_t *a_v;
                RkMatrix<HMAT_CHAMELEON_Z_t> *rkmat = node->rk();
                ScalarArray<HMAT_CHAMELEON_Z_t> *a_arr = rkmat->a;
                ScalarArray<HMAT_CHAMELEON_Z_t> *b_arr = rkmat->b;
                const rpk_ctx_t *ctx = runtime_rpk_zctx_get();
                rpk_matrix_t *Bra;
                assert( rank > 0 );
                assert( sizeof(HMAT_CHAMELEON_Z_t) == sizeof(CHAMELEON_Complex64_t) );
                a_u = (CHAMELEON_Complex64_t*)a_arr->ptr();
                a_v = (CHAMELEON_Complex64_t*)b_arr->ptr();
            
                Bra = (rpk_matrix_t*)CHAM_tile_get_ptr( tileB );
                
                /* Resize the rapack matrix to accomodate the new rank */
                rpkx_zlrsze(
                    ctx,
                    0,
                    tileB->m, tileB->n,
                    Bra,
                    rank, rank,
                    rank
                );

                assert( Bra->rk == rank);
                assert( Bra->rkmax == rank);
                assert( Bra->u );
                assert( Bra->v );

                /* Set the content of the new rapack matrix */
                CORE_zlacpy(
                    ChamUpperLower,
                    tileA->m, rank,
                    a_u, a_arr->lda,
                    (CHAMELEON_Complex64_t*)Bra->u, tileB->m
                );

                /* In the Hmat-oss library, v is stored non-transposed
                 * so we need to transpose it before copying it to the rapack matrix */
                CORE_zgeadd(
                    ChamTrans,
                    rank, tileA->n,
                    1.0,
                    a_v, b_arr->lda,
                    0.0,
                    (CHAMELEON_Complex64_t*)Bra->v, Bra->rkmax
                );
            }
        }
        
        virtual ~CHAMELEON_ZHmat2Rpk_Proc() {

        }
    };
}

#endif

void
TCORE_zhmat2rpk( CHAM_tile_t *tileA, CHAM_tile_t *tileB )
{

    assert( tileA->m == tileB->m );
    assert( tileA->n == tileB->n );

    hmat_matrix_t *Ah = (hmat_matrix_t*)CHAM_tile_get_ptr( tileA );
    rpk_matrix_t *Bra = (rpk_matrix_t*)CHAM_tile_get_ptr( tileB );
    hmat_interface_t *hi = (hmat_interface_t*)runtime_hmat_zinterface_get();
    assert( hi );

    CHAMELEON_ZHmat2Rpk_Proc tree_proc(tileA, tileB);
    hmat_procedure_t hproc;

    /* Walk the hierarchical matrix to fill the rapack matrix */
    hproc.value_type = HMAT_DOUBLE_COMPLEX;
    hproc.internal = (void*)&tree_proc;
    hi->walk(Ah, &hproc);
}

