#ifndef LOCALITY_SENSITIVE_HASHING
#define LOCALITY_SENSITIVE_HASHING

#include<vector>
#include<unordered_set>
#include<unordered_map>

template< typename real >
void sample_normal( real &x, real &y ){
    real u=drand48(), v=drand48();
    real lnu = sqrt( -2.0*log(u) );
    x = lnu*cos(2.0*M_PI*v);
    y = lnu*sin(2.0*M_PI*v);
}

template< int NUM_TABLES=10, typename real=float, typename stored=int >
class locality_sensitive_hash {
private:
    std::vector< real >                     m_hash_vec[NUM_TABLES];
    //real                                    m_hash_off[NUM_TABLES];
    std::unordered_multimap< int, stored >  m_table[NUM_TABLES];
    real                                    m_r;
    int                                     m_nd;
    
    inline int hash_value_func( const int fid, const std::vector<real> &coord ) const {
        real sum = 0.0;
        for( int i=0; i<m_nd; i++ ){
            sum += m_hash_vec[fid][i]*coord[i];
        }
        sum += (drand48()-0.5)*m_r;
        return int(floor(sum/m_r));
    }

    void initialize(){
        real x,y;
        for( int i=0; i<NUM_TABLES; i++ ){
            //m_hash_off[i] = drand48()*m_r;
            m_hash_vec[i].resize( m_nd );
            real len = 0.0;
            for( int j=0; j<m_nd; j++ ){
                sample_normal( x, y );
                m_hash_vec[i][j] = x; //drand48() > 0.5 ? x : -x;
                len += x*x;
            }
            len = sqrt(len);
            for( int j=0; j<m_nd; j++ ){
                m_hash_vec[i][j] /= len;
            }
        }
    }

public:
    locality_sensitive_hash( const int nd, const real r ) {
        m_nd = nd;
        m_r = r;
        initialize();
    }

    void insert( const stored &id, const std::vector<real> &coord ){
        for( int i=0; i<NUM_TABLES; i++ ){
            int hv = hash_value_func( i, coord );
            m_table[i].insert( std::pair< int, stored >( hv, id ) );
        }
    }
    
    std::unordered_set<stored> query_int( const int fid, const std::vector<real> &coord ) const {
        std::unordered_set<stored> results;
        for( int i=0; i<20; i++ ){
            int hv = hash_value_func( fid, coord );
            auto range = m_table[fid].equal_range( hv );
            for( auto it=range.first; it!=range.second; ++it ){
                results.insert( it->second );
            }
        }
        return results;
    }

    std::vector< stored > query( const std::vector<real> &coord ) const {
        std::vector< stored > out;
        std::unordered_set< stored > results = query_int( 0, coord );
        for( int i=1; i<NUM_TABLES; i++ ){
            std::unordered_set< stored > tmp = query_int( i, coord );
            for( auto it=tmp.begin(); it!=tmp.end(); ++it ){
                results.insert( *it );
            }
        }
        for( auto it=results.begin(); it!=results.end(); ++it ){
            out.push_back( *it );
        }
        return out;
    }
    
    std::vector< stored> query2( const std::vector<real> &coord ) const {
        std::vector< stored > out;
        std::unordered_set< stored > exclude, results = query_int( 0, coord );
        for( int i=1; i<NUM_TABLES; i++ ){
            std::unordered_set< stored > tmp = query_int( i, coord );
            for( auto it=results.begin(); it!=results.end(); ++it ){
                if( tmp.find( *it ) == tmp.end() )
                    exclude.insert( *it );
            }
        }
        for( auto it=results.begin(); it!=results.end(); ++it ){
            if( exclude.find( *it ) == exclude.end() ){
                out.push_back( *it );
            }
        }
        return out;
    }

};

#endif