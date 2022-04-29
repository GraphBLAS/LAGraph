import re
import os
import sys

def usage():
    print("func_extract.py <HEADER_FILE_PATH>")

def main():
    if len(sys.argv) != 2:
        usage()
        exit(1)

    header_filename = sys.argv[1]

    # make sure input file exists
    if not os.path.isfile(header_filename):
        print("{} seems to be missing\n".format(header_filename))
        exit(1)

    # regex pattern
    # looking for:
    # "LAGRAPH_PUBLIC" followed by new line
    # followed by anything (can be multiple lines)
    # ends with new line ) ;
    #
    # e.g. 
    #
    # LAGRAPH_PUBLIC
    # int LAGraph_TriangleCount
    # (
    #    // output:
    #    uint64_t      *ntriangles,   // # of triangles
    #    // input/output:
    #    LAGraph_Graph  G,
    #    char          *msg
    # ) ;
    #
    pattern = "LAGRAPH_PUBLIC\n((.|\n)*?)\n\) ;"
    pattern_enum = "typedef enum\n((.|\n)*?) ;"

    with open(header_filename, "r") as f:
        header = f.read()
        matches = re.findall(pattern, header)
        for func in matches:
            print("{0}\n\n".format("\\begin{verbatim}\n" + func[0] +
            "\n) ;\n\\end{verbatim}\n\n"))
        matches = re.findall(pattern_enum, header)
        for func in matches:
            print("{0}\n\n".format("\\begin{verbatim}\ntypedef enum\n"
            + func[0] + " ;\n\\end{verbatim}\n\n"))

if __name__ == "__main__":
    main()

