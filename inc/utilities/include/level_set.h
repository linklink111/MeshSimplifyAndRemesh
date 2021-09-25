#ifndef LEVEL_SET_H
#define LEVEL_SET_H

#include<algorithm>

template< typename real >
real ls_abs( const real &in ){
    return (in < real(0.0) ? -in : in);
}

template< typename real, typename real3, typename distance_func >
real3 ls_signed_distance_gradient( distance_func &f, const real3 &query, const real eps=1e-3 ){
    const real f0 = f( query );
    const real f1 = f( query+real3(eps,0,0) );
    const real f2 = f( query+real3(0,eps,0) );
    const real f3 = f( query+real3(0,0,eps) );
    return real3( f1-f0, f2-f0, f3-f0 )/eps;
}

template< typename real, typename real3 >
real ls_signed_distance_point( const real &point, const real3 &query ){
    return query-point.length();
}

// negative inside the sphere
template< typename real, typename real3 >
real ls_signed_distance_sphere( const real3 &center, const real r, const real3 &query ){
    return (query-center).length()-r;
}

// negative inside the box
template< typename real, typename real3 >
real ls_signed_distance_box( const real3 &minim, const real3 &maxim, const real3 &query ){
    real3 clamp = query;
    if( query[0] >= minim[0] && query[0] <= maxim[0] && query[1] >= minim[1] && query[1] <= maxim[1] && query[2] >= minim[2] && query[2] <= maxim[2] ){
        // inside the box
        clamp[0] = std::min( ls_abs(query[0]-minim[0]), ls_abs(query[0]-maxim[0]) );
        clamp[1] = std::min( ls_abs(query[1]-minim[1]), ls_abs(query[1]-maxim[1]) );
        clamp[2] = std::min( ls_abs(query[2]-minim[2]), ls_abs(query[2]-maxim[2]) );
        return - std::min( clamp[0], std::min( clamp[1], clamp[2] ) );
    }
    // outside the box, can clamp the query point to the
    // bounding box to get the closest point (since the
    // box is orthogonal).
    clamp[0] = std::max( minim[0], std::min( maxim[0], query[0] ) );
    clamp[1] = std::max( minim[1], std::min( maxim[1], query[1] ) );
    clamp[2] = std::max( minim[2], std::min( maxim[2], query[2] ) );
    return (query-clamp).length();
}



#endif