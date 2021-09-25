#ifndef OPTIONS_H
#define OPTIONS_H

#include<map>
#include<string>

#include<string_utils.h>

namespace utilities {

    class options {
    private:
        std::map< std::string, std::string > m_options;
    public:
        
        std::string &operator[]( const std::string &in ){
            return m_options[in];
        }
        
        bool exists( const std::string &name ){
            return m_options.find( name ) != m_options.end();
        }
        
        template< typename T >
        void set_default( const std::string &name, const T &value ){
            if( !exists(name) )
                m_options[name] = to_str(value);
        }
        
        template< typename T >
        void set( const std::string &name, const T &value ){
            m_options[name] = to_str(value);
        }
        
        template< typename T >
        T get( const std::string &name ){
            return from_str<T>( m_options[name] );
        }
    };
    
};

#endif