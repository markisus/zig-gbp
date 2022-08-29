import re

api_path = "/home/mark/blasfeo/include/blasfeo_d_blasfeo_api.h"
pattern = r"^(double|void) blasfeo_(?P<funcname>[a-zA-Z_]+)[(](?P<args>[^)]*)[)]"

def try_parse_basic_index(argname):
    if argname[-1] in ["i", "j"]:
        # print("parsed index", argname)
        return {
            'refers_to': argname[:-1],
            'ij': argname[-1]
        }

def try_parse_index(argname):
    # ending with i or j
    basic_parse = try_parse_basic_index(argname)
    if basic_parse:
        return basic_parse
    
    if argname[-1] in ["0", "1", "2", "3"]:
        # print("have index", argname)
        potential_parse = try_parse_basic_index(argname[:-1])
        if not potential_parse:
            return None
        potential_parse['refers_to'] += argname[-1]
        # print("returning", potential_parse)
        return potential_parse
        
    elif len(argname) > 2 and argname[-2:] in ["_t", "_n"]:
        potential_parse = try_parse_basic_index(argname[:-2])
        if not potential_parse:
            return None
        potential_parse['refers_to'] += argname[-2:]
        return potential_parse

def analyze_declaration(args):
    # search arguments for matrices
    arginfos = {}
    matrix_args = []
    vector_args = []

    # first get all the matrix and vector inputs
    for argtype, argname in args:
        if argtype == "blasfeo_dmat":
            # argname is like "*sA"
            assert argname[:2] == "*s"
            arginfos[argname] = {
                'type': 'matrix',
                'zig': argname[2:]
            }
            matrix_args.append(argname[2:])
        elif argtype == "blasfeo_dvec":
            # argname is like "*sx"
            assert argname[:2] == "*s"
            arginfos[argname] = {
                'type': 'vector',
                'zig': argname[2:]
            }
            vector_args.append(argname[2:])

    # then find all the indexing / dimension variables
    for argtype, argname in args:
        if argtype == "int":
            if argname == "m" or argname == "kmax":
                if matrix_args:
                    matrix_name = matrix_args[0]
                    arginfos[argname] = {
                        'type': 'dim',
                        'zig': f'{matrix_name}.rows'
                    }
                else:
                    vector_name = vector_args[0]
                    arginfos[argname] = {
                        'type': 'dim',
                        'zig': f'{vector_name}.len'
                    }
            elif argname == "n":
                matrix_name = matrix_args[0]
                arginfos[argname] = {
                    'type': 'dim',
                    'zig': f'{matrix_name}.cols'
                }
            elif argname == "k":
                matrix_name = matrix_args[-1]
                arginfos[argname] = {
                    'type': 'dim',
                    'zig': f'{matrix_name}.cols'
                }
            else:
                index_info = try_parse_index(argname)
                if not index_info:
                    continue
                if index_info['ij'] == 'i':
                    name = index_info['refers_to']
                    if name.upper() in matrix_args:
                        arginfos[argname] = {
                            'type': 'offset',
                            'zig': f'{name.upper()}.row_start'
                        }
                    elif name in vector_args:
                        arginfos[argname] = {
                            'type': 'offset',
                            'zig': f'{name}.start'
                        }
                else:
                    assert index_info['ij'] == 'j'
                    name = index_info['refers_to']
                    if name.upper() in matrix_args:                    
                        arginfos[argname] = {
                            'type': 'offset',
                            'zig': f'{name.upper()}.col_start'
                        }
    return arginfos

def translate_declaration(funcname, returntype, args, arginfos):
    paramlist = []
    for argtype, argname in args:
        arginfo = arginfos.get(argname, None)
        if not arginfo:
            zigtype = 'c_int'
            if argtype == 'double':
                zigtype = 'f64'
            elif argtype == 'void':
                zigtype = 'anyopaque'
            if argname.startswith("*"):
                argname = argname[1:]
                zigtype = "*" + zigtype
            paramlist.append(f'{argname} : {zigtype}')
        elif arginfo['type'] == 'dim':
            # skip dimension variables
            continue
        elif arginfo['type'] == 'matrix':
            paramlist.append(f'{arginfo["zig"]} : Matrix')
        elif arginfo['type'] == 'vector':
            paramlist.append(f'{arginfo["zig"]} : Vector')
        else:
            assert arginfo['type'] == 'offset'
            # paramlist.append(f'{argname} : c_int')

    params = ", ".join(paramlist)
    return f'pub fn {funcname}({params}) {returntype}'

def get_wrapped_impl(funcname, args, arginfos):
    zigvars_list = []
    for argtype, argname in args:
        arginfo = arginfos.get(argname, None)
        if arginfo:
            if arginfo['type'] == 'matrix':
                zigvars_list.append(f'{arginfo["zig"]}.dmat')
            elif arginfo['type'] == 'vector':
                zigvars_list.append(f'{arginfo["zig"]}.dvec')
            else:
                zigvars_list.append(f'{arginfo["zig"]}')
        else:
            if argname[0] == "*":
                argname = argname[1:]
            zigvars_list.append(argname)

    zigvars = ", ".join(zigvars_list)
    return f'{{ cblasfeo.blasfeo_d{funcname}({zigvars}); }}'

with open(api_path, "r") as f:
    for line in f.readlines():
        match = re.match(pattern, line)

        if not match:
            continue

        returntype = match.group(1)
        funcname = (match.group('funcname'))
        if funcname.startswith('cm_'):
            # skip complex
            continue

        assert funcname.startswith('d')
        funcname = funcname[1:]
        
        args = [a.strip() for a in (match.group('args')).split(",")]
        args = [a.split(" ")[-2:] for a in args]

        arginfos = analyze_declaration(args)
        translated_declaration = translate_declaration(
            funcname, returntype, args, arginfos)

        print(translated_declaration)
        print(get_wrapped_impl(funcname, args, arginfos))

        # # search arguments for matrices
        # matrix_args = []
        # vector_args = []
        # for argtype, argname in args:
        #     if argtype == "blasfeo_dmat":
        #         # argname is like "*sA"
        #         assert argname[:2] == "*s"
        #         matrix_args.append(argname[2:].lower())
        #     if argtype == "blasfeo_dvec":
        #         # argname is like "*sx"
        #         assert argname[:2] == "*s"
        #         vector_args.append(argname[2:].lower())

        # output_args = []
        # for argtype, argname in args:
        #     if argname.startswith('*s'):
        #         argname = argname[2:]
        #     if argtype == 'int':
        #         # this arg is related to a dimension, skip
        #         continue
        #     elif argtype == 'blasfeo_dvec':
        #         output_args.append(f"{argname.lower()} : Vector")
        #     elif argtype == 'blasfeo_dmat':
        #         variable = argname[-1]
        #         output_args.append(f"{argname.upper()} : Matrix")
        #     elif argtype == 'double':
        #         output_args.append(f"{argname} : f64")

        # output_returntype = "f64" if returntype == "double" else "void"
        # output_decl = f"pub fn {funcname}({', '.join(output_args)}) {output_returntype}"

        # # function body:
        # zigvars = []
        # for argtype, argname in args:
        #     zigvar = argname
        #     if zigvar.startswith('*s'):
        #         zigvar = argname[2:]
        #     if argtype == 'blasfeo_dvec':
        #         zigvar = zigvar.lower() + '.dvec'
        #     elif argtype == 'blasfeo_dmat':
        #         zigvar = zigvar.upper() + '.dmat'
        #     elif argtype == 'int':
        #         # this arg is related to a dimension, skip
        #         if zigvar == 'm':
        #             if matrix_args:
        #                 zigvar = matrix_args[0] + '.rows'
        #             else:
        #                 zigvar = vector_args[0] + '.len'
        #         elif zigvar == 'n':
        #             zigvar = matrix_args[0].upper() + '.cols'
        #         elif zigvar == 'k':
        #             zigvar = matrix_args[-1].upper() + '.cols'
        #         elif zigvar.endswith('i'):
        #             zigvar = zigvar[:-1]
        #             if zigvar in matrix_args:
        #                 zigvar = zigvar.upper() + '.row_start'
        #             else:
        #                 assert zigvar in vector_args
        #                 zigvar = zigvar + '.start'
        #         elif zigvar.endswith('j'):
        #             zigvar = zigvar[:-1]
        #             assert zigvar in matrix_args
        #             zigvar = zigvar.upper() + '.col_start'
        #         else:
        #             # print("unhandled int", argname)
        #             continue
        #     elif argtype == 'double':
        #         pass
        #     zigvars.append(zigvar)                
            
        
        # wrappedcall = f"cblasfeo.blasfeo_d{funcname}({', '.join(zigvars)})"

        # print(f"{output_decl} {{\n\t{wrappedcall};\n}}")
        # # call the cimpl
        

        # # parse line
        # # void basfeo_[funcname](params)
        # # int m => A.rows
        # # int n => A.cols
        # # basfeo_dmat *sA => A.dmat
        # # int ai => A.row_start
        # # int aj => A.col_start
        # # basfeo_dvec *sx => x.dvec
        # # int xi => x.start
        # # blasfeo_dvec *sz => z.dvec
        # # int zi => zi.start

        # # 
        # # print(line)
