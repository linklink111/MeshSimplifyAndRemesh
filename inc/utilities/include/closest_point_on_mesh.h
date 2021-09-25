#ifndef CLOSEST_POINT_ON_MESH_H
#define CLOSEST_POINT_ON_MESH_H

#include<vector>
#include<algorithm>

#include"octree.h"
#include"triangle.h"

namespace geom {
    
    
    template< typename real, typename real3 >
    class closest_point_on_mesh {
    private:
        real3                   m_minim;
        real3                   m_maxim;
        std::vector< real3 >    m_verts;
        std::vector< real3 >    m_normal;
        std::vector< int>       m_tris;
        geom::octree_node<real,real3,int>      m_octree;
        
        inline real3 safe_normalize( const real3 &v, const real eps=real(1e-10) ) const {
            real L = std::max( eps, v.length() );
            return v/L;
        }
        
        inline real angle( const real3 &A, const real3 &B, const real3 &C ) const {
            real3 u = safe_normalize( A-B );
            real3 v = safe_normalize( C-A );
            return acos( u.dot(v) );
        }
    public:
        closest_point_on_mesh(){
            
        }
        
        void initialize( const std::vector<real3> & verts, const std::vector<int> & tris, int subdiv_level=5 ){
            m_verts = verts;
            m_tris  = tris;
            if( m_verts.size() == 0 ){
                m_minim = real3(0,0,0);
                m_maxim = real3(1,1,1);
                return;
            }
            m_normal.resize( m_verts.size() );
            m_minim = m_verts[0];
            m_maxim = m_verts[0];
            m_normal[0] = real3(0,0,0);
            for( int i=1; i<m_verts.size(); i++ ){
                m_minim = m_minim.min( m_verts[i] );
                m_maxim = m_maxim.max( m_verts[i] );
                m_normal[i] = real3(0,0,0);
            }
            m_octree = geom::octree_node<real,real3,int>(m_minim,m_maxim);
            m_octree.subdivide_n_levels( subdiv_level );
            for( int i=0; i<m_tris.size(); i+=3 ){
                real3 A = m_verts[m_tris[i+0]];
                real3 B = m_verts[m_tris[i+1]];
                real3 C = m_verts[m_tris[i+2]];
                real3 N = safe_normalize( (B-A).cross(C-A) );
                m_normal[m_tris[i+0]] += N*angle( C, A, B );
                m_normal[m_tris[i+1]] += N*angle( A, B, C );
                m_normal[m_tris[i+2]] += N*angle( B, C, A );
                real3 triminim = A;
                real3 trimaxim = A;
                triminim = triminim.min( B );
                trimaxim = trimaxim.max( B );
                triminim = triminim.min( C );
                trimaxim = trimaxim.max( C );
                //m_octree.add_item( i/3, triminim, trimaxim );
                try {
                    m_octree.add_item_predicate( *this, i/3 );
                } catch( const char *e ){}
            }
        }
        
        bool operator()( int id, const real3& minim, const real3& maxim ){
            int vid[] = { m_tris[id*3+0], m_tris[id*3+1], m_tris[id*3+2] };
            real3 cen = (maxim+minim)/2.0;
            real3 n, cp = triangle_closest_point<real,real3>( m_verts[vid[0]], m_normal[vid[0]], m_verts[vid[1]], m_normal[vid[1]], m_verts[vid[2]], m_normal[vid[2]], cen, n );
            //return cp[0] >= minim[0] && cp[0] <= maxim[0] && cp[1] >= minim[1] && cp[1] <= maxim[1] && cp[2] >= minim[2] && cp[2] <= maxim[2];
            return (cp-cen).length_squared() < ((maxim-minim)/2.0).length_squared();
        }
        
        real3 operator()( const real3 &p, const real size, real3 &n ) const {
            std::vector<int> tris;
            real3 minim( p[0]-size, p[1]-size, p[2]-size );
            real3 maxim( p[0]+size, p[1]+size, p[2]+size );
            m_octree.get_items_overlapping_aabb( minim, maxim, tris );
            
            real td, mind = 1e10;
            real3 closest_point, closest_normal, tmp, tnormal;
            for( int i=0; i<tris.size(); i++ ){
                int vid[] = { m_tris[tris[i]*3+0], m_tris[tris[i]*3+1], m_tris[tris[i]*3+2] };
                try {
                    tmp = triangle_closest_point<real,real3>( m_verts[vid[0]], m_normal[vid[0]], m_verts[vid[1]], m_normal[vid[1]], m_verts[vid[2]], m_normal[vid[2]], p, tnormal );
                td = (p-tmp).length_squared();
                if( td < mind ){
                    mind = td;
                    closest_point = tmp;
                    closest_normal = tnormal;
                }
                } catch( const char* ){}
            }
            n = closest_normal;
            return closest_point;
        }
    
    };
    
};

#endif