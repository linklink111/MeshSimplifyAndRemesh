#ifndef COMMAND_LINE_ARGUMENTS_H
#define COMMAND_LINE_ARGUMENTS_H

#include<ctype.h>
#include<cstdio>
#include<cstdarg>
#include<cstdlib>
#include<cassert>

/**
 @file command_line_arguments.h
 @author James Gregson james.gregson@gmail.com
 @brief provides some basic command line argument handling.  Free for all use. Sample usage below:
 
 #include"command_line_arguments.h"
 
 char input_file[1024];
 char output_file[1024];
 int  index_offset = 0;
 bool output_tris = false;
 bool output_quads = false;
 bool output_tets = true;
 bool output_hexes = true;
 bool output_surface = true;
 
 bool load_arguments( int argc, char **argv ){
    bool ok = true;
    ok &= command_line_get_argument( argc, argv, "-input", ARGUMENT_STRING, input_file );
    ok &= command_line_get_argument( argc, argv, "-output", ARGUMENT_STRING, output_file );
    command_line_get_argument( argc, argv, "-index_offset", ARGUMENT_INT, &index_offset );
    command_line_get_argument( argc, argv, "-save_tris", ARGUMENT_BOOL, &output_tris );
    command_line_get_argument( argc, argv, "-save_quads", ARGUMENT_BOOL, &output_quads );
    command_line_get_argument( argc, argv, "-save_tets", ARGUMENT_BOOL, &output_tets );
    command_line_get_argument( argc, argv, "-save_hexes", ARGUMENT_BOOL, &output_hexes );
    command_line_get_argument( argc, argv, "-save_surface", ARGUMENT_BOOL, &output_surface );
    if( !ok ){
        std::cout << "Convert .mesh to .vtu" << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "  " << argv[0] << " [REQUIRED] [OPTIONAL]" << std::endl;
        std::cout << "[REQUIRED] - All of:" << std::endl;
        std::cout << "  -input        <filename>   input filename" << std::endl;
        std::cout << "  -output       <filename>   output filename" << std::endl;
        std::cout << "[OPTIONAL] - Any combination of:" << std::endl;
        std::cout << "  -index_offset <integer>    positive or negative offset added to nodal ID when writing output " << std::endl;
        std::cout << "  -save_tris    <true|false> save any triangles found in input" << std::endl;
        std::cout << "  -save_quads   <true|false> save any triangles found in input" << std::endl;
        std::cout << "  -save_tets    <true|false> save any tetrahedra found in input" << std::endl;
        std::cout << "  -save_hexes   <true|false> save any hexahedra found in input" << std::endl;
        std::cout << "  -save_surface <true|false> save all surface elements" << std::endl;
        exit(1);
    }
    return ok;
 }

*/

/**
 Defines the possible types for a command line argument
*/
typedef enum {
    ARGUMENT_FLAG,
	ARGUMENT_BOOL,
    ARGUMENT_INT,
    ARGUMENT_FLOAT,
	ARGUMENT_DOUBLE,
    ARGUMENT_STRING,
    ARGUMENT_INPUT_FILE,
} arguments_type;

/**
 Error function for command line parsing error
*/
static void command_line_error( const char *file, int line, const char *fmt, ... ){
    va_list args;
    va_start (args, fmt);
    fprintf( stderr, "Error at %s line %d: ", file, line );
    vfprintf( stderr, fmt, args );
    va_end (args );
    assert( false );
}

/**
 Get a command-line argument by searching through the argument list.  Ensures that the 
 argument was properly converted by sscanf, and for ARGUMENT_INPUT_FILE, checks that 
 the file can be opened for reading.
*/
static bool command_line_get_argument( const int argc, const char **argv, const char *arg_name, arguments_type arg_type, void *arg_data, int count=1 ){
    FILE *fp=NULL;
    int k;
    for( int i=1; i<argc; i++ ){
        if( strcmp( arg_name, argv[i] ) == 0 && (arg_type == ARGUMENT_FLAG || i+count < argc) ){
            switch( arg_type ){
                case ARGUMENT_FLAG:
                    *((bool*)arg_data) = true;
                    break;
				case ARGUMENT_BOOL:
                    for( k=0; k<count; k++ ){
                        if( strcmp( argv[i+1+k], "true" ) == 0 ){
                            *((bool*)arg_data) = true;
                            return true;
                        } else if( strcmp( argv[i+1+k], "false" ) == 0 ){
                            *((bool*)arg_data) = false;
                            return true;
                        }
                        command_line_error( __FILE__, __LINE__, "could not parse boolean argument for argument %s from %s\n", arg_name, argv[i+1+k] );
                    }
					return false;
					break;
                case ARGUMENT_INT:
                    for( k=0; k<count; k++ ){
                        if( sscanf( argv[i+1+k], "%d", (int*)arg_data+k ) != 1 ){
                            command_line_error( __FILE__, __LINE__, "could not parse integer argument for argument %s from %s\n", arg_name, argv[i+1+k] );
                            return false;
                        }
                    }
                    return true;
                    break;
                case ARGUMENT_FLOAT:
                    for( k=0; k<count; k++ ){
                        if( sscanf( argv[i+1+k], "%f", (float*)arg_data+k ) != 1 ){
                            command_line_error( __FILE__, __LINE__, "could not parse floating point argument for argument %s from %s\n", arg_name, argv[i+1+k] );
                            return false;
                        }
                    }
                    return true;
                    break;
                case ARGUMENT_DOUBLE:
                    for( k=0; k<count; k++ ){
                        if( sscanf( argv[i+1+k], "%lf", (double*)arg_data+k ) != 1 ){
                            command_line_error( __FILE__, __LINE__, "could not parse double precision argument for argument %s from %s\n", arg_name, argv[i+1+k] );
                            return false;
                        }
                    }
                    return true;
                    break;
                case ARGUMENT_STRING:
                    if( sscanf( argv[i+1], "%s", (char*)arg_data ) != 1 ){
                        command_line_error( __FILE__, __LINE__, "could not parse string argument for argument %s from %s\n", arg_name, argv[i+1] );
                        return false;
                    }
                    return true;
                    break;
                case ARGUMENT_INPUT_FILE:
                    if( sscanf( argv[i+1], "%s", (char*)arg_data ) != 1 ){
                        command_line_error( __FILE__, __LINE__, "could not parse input filename for argument %s from %s\n", arg_name, argv[i+1] );
                        return false;
                    }
                    fp = fopen( (char*)arg_data, "r" );
                    if( !fp ){
                        command_line_error( __FILE__, __LINE__, "input file %s for argument %s does not exist (or cannot be opened for reading\n", argv[i+1], arg_name );
                        return false;
                    }
                    fclose(fp);
                    return true;
                    break;
            }
        }
    }
    // argument was not found
    return false;
}

#endif