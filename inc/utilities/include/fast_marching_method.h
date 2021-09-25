#ifndef FAST_MARCHING_METHOD_H
#define FAST_MARCHING_METHOD_H

/**
 @module    fast_marching_method
 @file      fast_marching_method.h
 @author    James Gregson (james.gregson@gmail.com)
 @copyright Use freely (but at your own risk) for commercial and non-commercial purposes provided this header remains intact.
 @brief     Implements a simple version of the Fast-Marching Method of Sethian and Osher for the construction of signed distance functions in 2D and 3D.  
*/

#include<cmath>
#include<set>
#include<algorithm>

/**
 @brief enumerated type to describe the state of a vertext in the domain
*/
typedef enum {
    FMM_FROZEN,             /** Vertex is frozen and has had its final value set */
    FMM_BOUNDARY_CONDITION, /** Vertex is an input boundary condition (special case of frozen) */
    FMM_ACTIVE,             /** Vertex is in the set of active vertices, but needs to have its velocity updated */
    FMM_FAR,                /** Vertex has not yet been reached by the algorithm (i.e. not in queue) */
    FMM_INVALID,
} fmm_state;

class fmm_callbacks {
public:
    virtual float       max_distance(){ return 1e10f; };
    virtual float       get_distance( int i, int j, int k ) = 0;
    virtual void        set_distance( int i, int j, int k, float val )=0;
    virtual fmm_state   get_state( int i, int j, int k )=0;
    virtual void        set_state( int i, int j, int k, fmm_state state )=0;
};

class fmm_pnt3i {
public:
    int x, y, z;
    fmm_pnt3i( int _x=0, int _y=0, int _z=0 ) : x(_x), y(_y), z(_z) {}
    
    fmm_pnt3i operator+( const fmm_pnt3i &in ) const {
        return fmm_pnt3i( x+in.x, y+in.y, z+in.z );
    }
    
    fmm_pnt3i operator-( const fmm_pnt3i &in ) const {
        return fmm_pnt3i( x-in.x, y-in.y, z-in.z );
    }
    
    bool operator<( const fmm_pnt3i &in ) const {
        return (x < in.x || (x==in.x && y<in.y) || (x==in.x && y==in.y && z<in.z));
    }
};

template< typename callbacks >
static void fast_marching_method_gradient_contribution( callbacks &cb, fmm_pnt3i &pnt, fmm_pnt3i &dir, float &A, float &B, float &C ){
    fmm_pnt3i    L=pnt-dir,                       R=pnt+dir;
    fmm_state    sL=cb.get_state(L.x,L.y,L.z),    sR=cb.get_state(R.x,R.y,R.z);
    float v, v2, dL=cb.get_distance(L.x,L.y,L.z), dR=cb.get_distance(R.x,R.y,R.z);
    if( sL == FMM_FROZEN && dL <= dR ){
        v = dL;
        L = L-dir;
        sL = cb.get_state(L.x,L.y,L.z);
        v2 = cb.get_distance(L.x,L.y,L.z);
        if( false && sL == FMM_FROZEN && v2 < v ){
            float a=+3.0/2.0;
            float b=-4.0/2.0;
            float c=+1.0/2.0;
            
            A += a*a;
            B += 2.0*( a*b*v + a*c*v2 );
            C += b*b*v*v + c*c*v2*v2 + 2.0*b*c*v*v2;
            
        } else {
            A += 1.0f;
            B -= 2.0f*v;
            C += v*v;
        }
        return;
    }
    if( sR == FMM_FROZEN && dR < dL ){
        v = dR;
        R = R+dir;
        sR = cb.get_state(R.x,R.y,R.z);
        v2 = cb.get_distance(R.x,R.y,R.z);
        if( false && sR == FMM_FROZEN && v2 < v ){
            float a=+3.0/2.0;
            float b=-4.0/2.0;
            float c=+1.0/2.0;
            
            A += a*a;
            B += 2.0*( a*b*v + a*c*v2 );
            C += b*b*v*v + c*c*v2*v2 + 2.0*b*c*v*v2;
        } else {
            A += 1.0f;
            B -= 2.0f*v;
            C += v*v;
        }
    }
    return;
}

template< typename callbacks >
float fast_marching_method_recompute_distance( callbacks &cb, fmm_pnt3i &pnt ){
    fmm_pnt3i dir[] = { fmm_pnt3i(1,0,0), fmm_pnt3i(0,1,0), fmm_pnt3i(0,0,1) };
    float A=0.0f, B=0.0f, C=-1.0f, det;
    fast_marching_method_gradient_contribution( cb, pnt, dir[0], A, B, C );
    fast_marching_method_gradient_contribution( cb, pnt, dir[1], A, B, C );
#if defined(FMM_USE_3D)
    fast_marching_method_gradient_contribution( cb, pnt, dir[2], A, B, C );
#endif
    
    det = B*B - 4.0*A*C;
    if( det < 0.0f ){
        //return cb.max_distance();
        return cb.get_distance( pnt.x, pnt.y, pnt.z );
    }
    det = sqrt( det );
    return (-B + det)/(2.0*A);
}

template< typename callbacks >
void fast_marching_method( callbacks &cb, int nseeds, int *seeds ){
    typedef std::pair<float,fmm_pnt3i> queue_entry;
    std::set< queue_entry > queue;
    float dist, new_dist;
    fmm_state state;
    fmm_pnt3i nbr;
    int id, nnbr = 4;
    
    const fmm_pnt3i nbr_off[] = {
        fmm_pnt3i( -1,  0,  0 ),
        fmm_pnt3i(  1,  0,  0 ),
        fmm_pnt3i(  0, -1,  0 ),
        fmm_pnt3i(  0,  1,  0 ),
        fmm_pnt3i(  0,  0, -1 ),
        fmm_pnt3i(  0,  0,  1 )
    };
    
    // initial conditions, add the seed points to the queue
    for( int i=0; i<nseeds; i++ ){
        id = i*3;
        queue.insert( queue_entry( cb.get_distance(seeds[id+0],seeds[i*3+1],seeds[id+2]), fmm_pnt3i( seeds[id+0],seeds[id+1], seeds[id+2] ) ) );
        cb.set_state( seeds[id+0], seeds[id+1], seeds[id+2], FMM_BOUNDARY_CONDITION );
    }
    
    // until we run out of entries to process
    while( !queue.empty() ){
        // grab the first entry and remove it from
        // the queue of points to process
        fmm_pnt3i pnt = queue.begin()->second;
        queue.erase( *queue.begin() );
        
        // check the state of the point, if it is already
        // frozen or is invalid, continue on in the loop
        // without processing the point further
        state = cb.get_state( pnt.x, pnt.y, pnt.z );
        if( state == FMM_FROZEN || state == FMM_INVALID )
            continue;
        
        // otherwise mark the point as frozen to prevent
        // further updates to its state
        if( state != FMM_BOUNDARY_CONDITION ){
            cb.set_distance( pnt.x, pnt.y, pnt.z, fast_marching_method_recompute_distance( cb, pnt ) );
        }
        cb.set_state( pnt.x, pnt.y, pnt.z, FMM_FROZEN );
        
        // loop over the neighboring points
        for( int i=0; i<nnbr; i++ ){
            
            // grab the neighboring point coordinates
            nbr = pnt + nbr_off[i];
            
            // grab the neighboring point state and if it is
            // either in the active list or marked as far,
            // perform more processing
            state = cb.get_state( nbr.x, nbr.y, nbr.z );
            if( state == FMM_FAR || state == FMM_ACTIVE ){
                
                // store the current distance of the point from
                // the zero set, then recompute the distance
                // based upon the updated neighboring distances
                dist = cb.get_distance( nbr.x, nbr.y, nbr.z );
                new_dist = fast_marching_method_recompute_distance( cb, nbr );
                
                // if the new distance is less than the old
                // distance, the point needs to be 'updated'
                if( new_dist < dist ){
                    
                    // first set the distance for the point, then
                    // insert it into the queue and mark it as active.
                    // the point is reinserted into the queue with the
                    // updated distance rather than having the queue
                    // updated. this allows a regular ordered set to be used
                    // instead of a custom heap class.  the asymptotic
                    // O(NlogN) order of the algorithm is unchanged since
                    // each point has only a constant number of neighbors and
                    // can only be added to the queue once per neighbor,
                    // however the queue can be up to nnbrs times as large as
                    // in the original FMM
                    cb.set_distance( nbr.x, nbr.y, nbr.z, new_dist );
                    queue.insert( queue_entry(new_dist,nbr) );
                    cb.set_state( nbr.x, nbr.y, nbr.z, FMM_ACTIVE );
                }
            }
        }
    }
}

#endif
