#ifndef DUAL_CONTOUR_H
#define DUAL_CONTOUR_H

#include<Eigen/Dense>

template< typename real, typename real3 >
real3 dual_contour( const int num_inter, const real3 *pi, const real3 *ni, const real3 &def, const real eps=0.01 ){
    Eigen::Matrix< double, 3, 3 > A;
    Eigen::Matrix< double, 3, 1 > b;
    real3 pcen(0,0,0);
    
    for( int i=0; i<3; i++ ){
        b[i] = 0.0;
        for( int j=0; j<3; j++ ){
            A(i,j) = 0.0;
        }
    }
    
    for( int k=0; k<num_inter; k++ ){
        const real3 &n = ni[k];
        pcen += pi[k]/real(num_inter);
        
        for( int i=0; i<3; i++ ){
            for( int j=0; j<3; j++ ){
                A(i,j) += n[i]*n[j];
            }
            b[i] += n[i]*(n.dot(pi[k]));
        }
    }
        
    if( true || fabs(A.determinant()) < 1e-3 ){
        A(0,0) += eps; b[0] += eps*def[0];
        A(1,1) += eps; b[1] += eps*def[1];
        A(2,2) += eps; b[2] += eps*def[2];
    }
    b = A.inverse()*b;
    return real3( b[0], b[1], b[2] );
}


template< typename real, typename real3 >
real3 dual_contour_hex2( const real3 *node, const real *dist, const real3 *grad, real eps=0.01, real *dual_dis=NULL ){
    const int hex_edge[][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},   // bottom edges
        {4, 5}, {5, 6}, {6, 7}, {7, 4},   // top edges
        {0, 4}, {1, 5}, {2, 6}, {3, 7},   // side edges
    };

    int s, e;
    real w;
    real3 dual_pos;
    int   int_cnt=0;
    real3 int_pnt[12];
    real3 int_nrm[12];
    for( int i=0; i<12; i++ ){
        s = hex_edge[i][0];
        e = hex_edge[i][1];
        if( (dist[s] < 0.0 && dist[e] >= 0.0) || (dist[s] >= 0.0 && dist[e] < 0.0) ){
            w = -dist[s]/(dist[e]-dist[s]);
            int_pnt[int_cnt] = (1.0-w)*node[s] + w*node[e];
            int_nrm[int_cnt] = (1.0-w)*grad[s] + w*grad[e];
            int_cnt++;
        }
    }
    dual_pos = ( node[0]+node[1]+node[2]+node[3]+node[4]+node[5]+node[6]+node[7] )*0.125;
    if( dual_dis )
        *dual_dis = (dist[0]+dist[1]+dist[2]+dist[3]+dist[4]+dist[5]+dist[6]+dist[7])*0.125;
    if( int_cnt > 0 ){
        dual_pos = dual_contour<real,real3>( int_cnt, int_pnt, int_nrm, dual_pos, eps );
        if(dual_dis)
            *dual_dis = 0.0;
    }
    return dual_pos;
}

#endif