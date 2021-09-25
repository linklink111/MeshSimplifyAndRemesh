#ifndef TRIANGLE_H
#define TRIANGLE_H

#include<cassert>
#include<algorithm>

#include<vec3.h>
//#include<geometric_predicates.h>

namespace geom {
  /*
    typedef enum {
        DOES_NOT_INTERSECT,
        DOES_INTERSECT,
        MAY_INTERSECT,
    } intersection_type;    
    
    template< typename real, typename real3 >
    inline intersection_type triangle_intersects_segment( const real3 &A, const real3 &B, const real3 &C, const real &P, const real &Q, const real eps=1e-10 ){
        static bool init=true;
        if( init ){
            init = false;
            exactinit();
        }

        const double a[]={A[0],A[1],A[2]};
        const double b[]={B[0],B[1],B[2]};
        const double c[]={C[0],C[1],C[2]};
        const double p[]={P[0],P[1],P[2]};
        const double q[]={Q[0],Q[1],Q[2]};
        
        double abcp, abcq, abpq, bcpq, capq;

        // check that the segment straddles the triangle plane
        abcp = orient3d( p, a, b, c );
        abcq = orient3d( q, a, b, c );
        if( (abcp >= 0 && abcq >= 0) || (abcp <= 0 && abcq <= 0) )
            return DOES_NOT_INTERSECT;
        
        // compute the orientation of the three edges
        abpq = orient3d( a, b, p, q );
        bcpq = orient3d( b, c, p, q );
        capq = orient3d( c, a, p, q );
        
        // check if the segment does not intersect any of the primitives
        if( abpq < 0 || bcpq < 0 || capq < 0 )
            return DOES_NOT_INTERSECT;
        
        // check if any of the orientations is ambiguous, in which case
        // return that the segment may intersect the primitive
        if( (abpq == 0 && bcpq >= 0 && capq >= 0) ||
            (abpq >= 0 && bcpq == 0 && capq >= 0) ||
            (abpq >= 0 && bcpq >= 0 && capq == 0) )
            return MAY_INTERSECT;
        
        return DOES_INTERSECT;
    }
   */
        
    
    template< typename real, typename real3 >
    inline real triangle_area_and_normal( const real3 &A, const real3 &B, const real3 &C, real3 &normal, const real eps=1e-12 ){
        real area;
        normal = (B-A).cross(C-A);
        area = normal.length();
        normal /= std::max( eps, area );
        return area*real(0.5);
    }
    
    template< typename real, typename real3 >
    inline real triangle_side( const real3 &A, const real3 &B, const real3 &C, const real3 &p ){
        return (p-A).dot( (B-A).cross(C-A) );
    }
    
    template< typename real, typename real3 >
    inline real3 plane_closest_point( const real3 &N, const real3 &ref, const real3 &p ){
        real lambda = N.dot(ref-p)/N.dot(N);
        return p + lambda*N;
    }
    
    template< typename real, typename real3 >
    inline real3 segment_barycentric_coords( const real3 &A, const real3 &B, const real3 &p, const real tol=1e-9 ){
        real3 n;
        real lambda, ntn, nta, ntp, ntb;
        
        // compute the non-normalized edge vector
        n = B-A;
        
        // check for a degenerate segment
        // in which case just return the
        // segment starting point
        ntn = n.dot(n);
        if( fabs(ntn) < tol )
            return A;
        
        // compute the plane intersections
        // with the edge vector and segment
        // start & end points
        nta = n.dot(A);
        ntb = n.dot(B);
        ntp = n.dot(p);
        
        // compute the offset along the edge
        // vector for the point
        lambda = (ntp - nta)/ntn;
        
        assert( fabs(n.dot(p)-n.dot(A*(1.0-lambda)+B*lambda) ) < 1e-5 );
        
        
        // return the barycentric coordinates
        return real3( 1.0-lambda, lambda );
    }
    
    template< typename real, typename real3 >
    inline real3 segment_closest_point( const real3 &A, const real3 &B, const real3 &p, const real tol=1e-9 ){
        real3 bary = segment_barycentric_coords( A, B, p, tol );
        
        real lambda = bary[1];
        lambda = std::max( 0.0, std::min( 1.0, lambda ) );
                
        // return the offset point
        return A + (B-A)*lambda;
    }
    
    template< typename real, typename real3 >
    inline real3 segment_closest_point( const real3 &A, const real3 &nA, const real3 &B, const real3 &nB, const real3 &p, real3 &normal, const real tol=1e-9 ){
        real3 bary = segment_barycentric_coords( A, B, p, tol );
        
        real lambda = bary[1];
        lambda = std::max( 0.0, std::min( 1.0, lambda ) );
        
        // return the offset point
        normal = (1.0-lambda)*nA + lambda*nB;
        return A + (B-A)*lambda;
    }
    
    
    
    template< typename real, typename real3 >
    inline real3 triangle_barycentric_coords( const real3 &A, const real3 &B, const real3 &C, const real3 &p, const real tol=1e-9 ){
        real3 u = A-C;
        real3 v = B-C;
        real3 off = p-C;
        
        real alpha, beta, gamma, det, rhs[2];
        real M[2][2], iM[2][2];
        
        // compute the normal equations for the closest point
        M[0][0] = u.dot(u); M[0][1] = u.dot(v); rhs[0] = u.dot(off);
        M[1][0] = v.dot(u); M[1][1] = v.dot(v); rhs[1] = v.dot(off);
        
        // check the determinant and bail if the triangle
        // is degenerate
        det = M[0][0]*M[1][1]-M[0][1]*M[1][0];
        if( fabs(det) < tol )
            throw "triangle_closest_point() error: degenerate triangle!";
        
        // compute the inverse of the 2x2 matrix
        iM[0][0] =  M[1][1]/det; iM[0][1] = -M[0][1]/det;
        iM[1][0] = -M[1][0]/det; iM[1][1] =  M[0][0]/det;
        
        // multiply by the inverse and compute the third weight
        alpha = iM[0][0]*rhs[0] + iM[0][1]*rhs[1];
        beta  = iM[1][0]*rhs[0] + iM[1][1]*rhs[1];
        gamma = 1.0-alpha-beta;

        return real3( alpha, beta, gamma );
    }
    
    template< typename real, typename real3 >
    inline real3 triangle_closest_point( const real3 &A, const real3 &B, const real3 &C, const real3 &p, const real tol=1e-9 ){
        // get the parametric coordinates
        
        real3 n = (B-A).cross(C-A);
        real3 cp = plane_closest_point<real,real3>( n, A, p );
        
        if( (cp-A).cross(B-A).dot(n) > 0.0 ){
            return segment_closest_point<real,real3>( A, B, p, tol );
        } else if( (cp-B).cross(C-B).dot(n) > 0.0 ){
            return segment_closest_point<real,real3>( B, C, p, tol );
        } else if( (cp-C).cross(A-C).dot(n) > 0.0 ){
            return segment_closest_point<real,real3>( C, A, p, tol );
        } else {
            return cp;
        }
    }
    
    template< typename real, typename real3 >
    inline real3 triangle_closest_point( const real3 &A, const real3 &nA, const real3 &B, const real3 &nB, const real3 &C, const real3 &nC, const real3 &p, real3 &normal, const real tol=1e-9 ){
        // get the parametric coordinates
        real3 bary = triangle_barycentric_coords( A, B, C, p, tol );
        real alpha=bary[0], beta=bary[1], gamma=bary[2];
        real3 cp = A*alpha + B*beta + C*gamma;
        real3 n = (B-A).cross(C-A);
        
        // check all three weights in turn, if any is negative then
        // the point is on the wrong side of the opposite segment,
        // i.e. if alpha < 0, then the point lies on the wrong side
        // of segment B,C, so return the closest point on that segment
        // otherwise the point is in the triangle or on a bounding
        // segment, so return the interpolated value
        if( (cp-A).cross(B-A).dot(n) > 0.0 ){
            return segment_closest_point<real,real3>( A, nA, B, nB, p, normal, tol );
        } else if( (cp-B).cross(C-B).dot(n) > 0.0 ){
            return segment_closest_point<real,real3>( B, nB, C, nC, p, normal, tol );
        } else if( (cp-C).cross(A-C).dot(n) > 0.0 ){
            return segment_closest_point<real,real3>( C, nC, A, nA, p, normal, tol );
        } else {
            normal = alpha*nA + beta*nB + gamma*nC;
            return cp;
        }
    }
    
};

#endif