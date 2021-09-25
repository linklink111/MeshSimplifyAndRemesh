// MeshSimplifyAndRemesh.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "Simplify.h"
#include <stdio.h>
#include <time.h> 
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <string_utils.h>
#include <command_line_options.h>

#include <remesher.h>

void showHelp(const char* argv[]) {
    const char* cstr = (argv[0]);
    printf("Usage: %s <input> <output> <ratio> <agressiveness)\n", cstr);
    printf(" Input: name of existing OBJ format mesh\n");
    printf(" Output: name for decimated OBJ format mesh\n");
    printf(" Ratio: (default = 0.5) for example 0.2 will decimate 80%% of triangles\n");
    printf(" Agressiveness: (default = 7.0) faster or better decimation\n");
    printf("Examples :\n");
#if defined(_WIN64) || defined(_WIN32)
    printf("  %s c:\\dir\\in.obj c:\\dir\\out.obj 0.2\n", cstr);
#else
    printf("  %s ~/dir/in.obj ~/dir/out.obj 0.2\n", cstr);
#endif
} //showHelp()

std::vector<std::string> tokenize( std::string &in ){
    std::vector<std::string> tokens;
    std::istringstream iss(in);
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string> >(tokens));
    return tokens;
}

template< typename real >
bool load_obj( const char *filename, std::vector<real>& coords, std::vector<int>& tris ){
    
    std::ifstream in(filename);
    if( !in.is_open() )
        return false;
    
    std::string line;
    while( std::getline(in,line) ){
        if( line.empty() )
            continue;
        std::istringstream iss(line);
        std::string token;
        
        iss >> token;
        if( token == "v" ){
            real x, y, z;
            iss >> x >> y >> z;
            coords.push_back( x );
            coords.push_back( y );
            coords.push_back( z );
        } else if( token == "f" ){
            std::vector<std::string> tokens = tokenize(line);
            if( tokens.size() != 4 ){
                std::cout << "input not a triangle mesh!" << std::endl;
                return false;
            }
            for( int i=1; i<tokens.size(); i++ ){
                std::replace( tokens[i].begin(), tokens[i].end(), '/', ' ');
                int id = from_str<int>(tokens[i])-1;
                tris.push_back( id );
            }
        }
    }
    
    return true;
}

int main(int argc, const char* argv[]) {
    char input_file[1024];
    char output_file_remesh[1024];
    char output_file_simp[1024];
    double edge_len = 1.0/50.0;
    double feat = acos(0.707)*180.0/M_PI;
    int iters=10;
    bool need_remesh = 0;
    bool need_simplify = 0;
    float reduceFraction = 0.5;
    double agressiveness = 7.0;
    
    utilities::command_line_options clopts;
    clopts.add_required_parameter( "-input",  ARGUMENT_STRING, input_file,  1, "input file in .obj format" );
    clopts.add_required_parameter( "-outputremesh", ARGUMENT_STRING, output_file_remesh, 1, "output file in .obj format for remeshing" );
    clopts.add_required_parameter( "-outputsimp", ARGUMENT_STRING, output_file_simp, 1, "output file in .obj format for simplification" );
    clopts.add_optional_parameter( "-size",   ARGUMENT_DOUBLE, &edge_len,   1, "edge length as fraction of bounding box longest side" );
    clopts.add_optional_parameter( "-feat",   ARGUMENT_DOUBLE, &feat,       1, "feature threshold, degrees" );
    clopts.add_optional_parameter( "-iters",  ARGUMENT_INT,    &iters,      1, "number of remeshing iterations to perform");
    clopts.add_optional_parameter( "-remesh",  ARGUMENT_BOOL,    &need_remesh,      1, "need to do remeshing");
    clopts.add_optional_parameter( "-simplify",  ARGUMENT_BOOL,    &need_simplify,      1, "need to do mesh simplifying");
    clopts.add_optional_parameter( "-ratio",  ARGUMENT_FLOAT,    &reduceFraction,      1, "reduce fraction for simplifying, 0.5 default");
    clopts.add_optional_parameter( "-agsv",  ARGUMENT_DOUBLE,    &agressiveness,      1, "agresiveness for simplifying, 7.0 default");
    if( !clopts.parse( argc, argv ) ){
        return 1;
    }
    // mesh simplification
    if(need_simplify){
        Simplify::load_obj(input_file);
        if ((Simplify::triangles.size() < 3) || (Simplify::vertices.size() < 3))
            return EXIT_FAILURE;
        int target_count = Simplify::triangles.size() >> 1;

        if (reduceFraction > 1.0) reduceFraction = 1.0; //lossless only
        if (reduceFraction <= 0.0) {
            printf("Ratio must be BETWEEN zero and one.\n");
            return EXIT_FAILURE;
        }
        target_count = round((float)Simplify::triangles.size() * reduceFraction);
        
        clock_t start = clock();
        printf("Input: %zu vertices, %zu triangles (target %d)\n", Simplify::vertices.size(), Simplify::triangles.size(), target_count);
        int startSize = Simplify::triangles.size();
        Simplify::simplify_mesh(target_count, agressiveness, true);
        //Simplify::simplify_mesh_lossless( false);
        if (Simplify::triangles.size() >= startSize) {
            printf("Unable to reduce mesh.\n");
            return EXIT_FAILURE;
        }
        Simplify::write_obj(output_file_simp);
        printf("Output: %zu vertices, %zu triangles (%f reduction; %.4f sec)\n", Simplify::vertices.size(), Simplify::triangles.size()
        , (float)Simplify::triangles.size() / (float)startSize, ((float)(clock() - start)) / CLOCKS_PER_SEC);

    }

    if(need_remesh){
        // remesh
        std::vector<double> coords, rm_coords;
        std::vector<int> tris, rm_tris;
        if( !load_obj( input_file, coords, tris ) )
            return 2;

        // compute the longest edge of the
        // input bounding box
        vec3d minim( coords[0], coords[1], coords[2] );
        vec3d maxim = minim;
        for( int i=3; i<coords.size(); i+=3 ){
            vec3d p( coords[i+0], coords[i+1], coords[i+2] );
            minim = minim.min( p );
            maxim = maxim.max( p );
        }
        vec3d delta = maxim-minim;
        double L = std::max( delta[0], std::max( delta[1], delta[2] ) );
        edge_len *= L;

        remesher_options opts;
        opts["REMESHER_REFINE_FEATURES"]        = "TRUE";
        opts["REMESHER_COARSEN_FEATURES"]       = "TRUE";
        opts["REMESHER_FEATURE_THRESHOLD"]      = to_str(cos(feat*M_PI/180.0));
        opts["REMESHER_MIN_EDGE_LENGTH"]        = to_str(edge_len);
        opts["REMESHER_MAX_EDGE_LENGTH"]        = to_str(edge_len*2.0);
        opts["REMESHER_TARGET_EDGE_LENGTH"]     = to_str(edge_len);
        opts["REMESHER_RELATIVE_EDGE_ERROR"]    = to_str(-1.0);
        opts["REMESHER_ITERATIONS"]             = to_str(iters);

        try {
            remesh( opts, coords, tris, rm_coords, rm_tris );
        } catch( const char *err ){
            std::cout << err << std::endl;
            return 1;
        }
        save_obj( output_file_remesh, rm_coords, rm_tris );

    }

    return EXIT_SUCCESS;
}


// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
