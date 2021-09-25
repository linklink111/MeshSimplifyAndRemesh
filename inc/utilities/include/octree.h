#ifndef OCTREE_H
#define OCTREE_H

#include<vector>

namespace geom {
    template< typename real, typename real3, typename stored >
    class octree_node {
    private:
        octree_node         *m_child;
        real3                m_minim;
        real3                m_maxim;
        std::vector<stored>  m_items;
    public:
        octree_node(){
            m_child = NULL;
        }
        
        octree_node( const real3 &minim, const real3 &maxim ){
            m_minim = minim;
            m_maxim = maxim;
            m_child = NULL;
        }
        
        ~octree_node(){
            if( m_child )
                delete [] m_child;
        }
        
        void set_extents( const real3 &minim, const real3 &maxim ){
            m_minim = minim;
            m_maxim = maxim;
        }
        
        void subdivide(){
            if( m_child ) return;
            m_child = new octree_node[8];
            real3 m_center = (m_minim + m_maxim)/2.0;
            m_child[0].set_extents( real3( m_minim[0],  m_minim[1],  m_minim[2] ), real3( m_center[0], m_center[1], m_center[2] ) );
            m_child[1].set_extents( real3( m_center[0], m_minim[1],  m_minim[2] ), real3( m_maxim[0],  m_center[1], m_center[2] ) );
            m_child[2].set_extents( real3( m_center[0], m_center[1], m_minim[2] ), real3( m_maxim[0],  m_maxim[1],  m_center[2] ) );
            m_child[3].set_extents( real3( m_minim[0],  m_center[1], m_minim[2] ), real3( m_center[0], m_maxim[1],  m_center[2] ) );
            
            m_child[4].set_extents( real3( m_minim[0],  m_minim[1],  m_center[2] ), real3( m_center[0], m_center[1], m_maxim[2] ) );
            m_child[5].set_extents( real3( m_center[0], m_minim[1],  m_center[2] ), real3( m_maxim[0],  m_center[1], m_maxim[2] ) );
            m_child[6].set_extents( real3( m_center[0], m_center[1], m_center[2] ), real3( m_maxim[0],  m_maxim[1],  m_maxim[2] ) );
            m_child[7].set_extents( real3( m_minim[0],  m_center[1], m_center[2] ), real3( m_center[0], m_maxim[1],  m_maxim[2] ) );
        }
        
        void subdivide_n_levels( int n ){
            if( m_child || n <= 0 ) return;
            subdivide();
            for( int i=0; i<8; i++ ){
                m_child[i].subdivide_n_levels( n-1 );
            }
        }
        
        void subdivide_to_item_count( const int count, const int max_levels=8 ){
            if( m_child || max_levels <= 0 || m_items.size() <= count ) return;
            subdivide();
            for( int i=0; i<8; i++ ){
                m_child[i].subdivide_to_item_count( count, max_levels-1 );
            }
        }
        
        void get_items( const real3 &pnt, std::vector<stored> &items ) const {
            if( pnt[0] >= m_minim[0] && pnt[1] >= m_minim[1] && pnt[2] >= m_minim[2] && pnt[0] <= m_maxim[0] && pnt[1] <= m_maxim[1] && pnt[2] <= m_maxim[2] ){
                if( m_child ){
                    for( int i=0; i<8; i++ ){
                        m_child[i].get_items( pnt, items );
                    }
                } else {
                    for( int i=0; i<m_items.size(); i++ ){
                        items.push_back( m_items[i] );
                    }
                }
            }
        }
        
        void get_items_overlapping_aabb( const real3 &minim, const real3 &maxim, std::vector<stored> &items ) const {
            if( minim[0] <= m_maxim[0] && maxim[0] >= m_minim[0] && minim[1] <= m_maxim[1] && maxim[1] >= m_minim[1] && minim[2] <= m_maxim[2] && maxim[2] >= m_minim[2] ){
                if( m_child ){
                    for( int i=0; i<8; i++ ){
                        m_child[i].get_items_overlapping_aabb( minim, maxim, items );
                    }
                } else {
                    for( int i=0; i<m_items.size(); i++ ){
                        items.push_back( m_items[i] );
                    }
                }
            }
        }
        
        template< typename predicate >
        void add_item_predicate( predicate& p, const stored item ){
            if( p(item,m_minim,m_maxim) ){
                if( m_child ){
                    for( int i=0; i<8; i++ ){
                        m_child[i].add_item_predicate( p, item );
                    }
                } else {
                    m_items.push_back( item );
                }
            }
        }
        
        void add_item( const stored item, const real3 &minim, const real3 &maxim ){
            if( minim[0] <= m_maxim[0] && minim[1] <= m_maxim[1] && minim[2] <= m_maxim[2] && maxim[0] >= m_minim[0] && maxim[1] >= m_minim[1] && maxim[2] >= m_minim[2] ){
                if( m_child ){
                    for( int i=0; i<8; i++ ){
                        m_child[i].add_item( item, minim, maxim );
                    }
                } else {
                    m_items.push_back( item );
                }
            }
        }
    };
    
};

#endif
