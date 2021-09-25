#ifndef SCOPE_H
#define SCOPE_H

#include<iostream>

namespace utilities {

#define utilities_scope_defines int utilities::scope::level=0;
  
    class scope {
    private:
		static int level;
    public:
		scope(){
			level++;
		}
		
		~scope(){
			level--;
		}
		
		inline int get_level() const {
			return level;
		}
    };
	
	std::ostream& operator<<( std::ostream &os, const scope &s ){
		for( int i=0; i<s.get_level(); i++ )
			os << "  ";
		return os;
	}
    
};


#endif