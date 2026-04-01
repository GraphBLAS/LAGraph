import argparse
import io
import re

import pcpp


class _Preprocessor(pcpp.Preprocessor):
    """Preprocessor that silently passes through unfound includes."""

    def on_include_not_found(
        self, is_malformed, is_system_include, curdir, includepath
    ):
        raise pcpp.OutputDirective(pcpp.Action.IgnoreAndPassThrough)

    def on_error(self, file, line, msg):
        pass  # suppress errors from missing include contents


def _preprocess(input_file):
    pp = _Preprocessor()
    pp.define("LG_JIT_KERNEL(name) __LG_JIT_KERNEL__(name)")
    pp.parse(open(input_file).read(), input_file)
    out = io.StringIO()
    pp.write(out)
    return out.getvalue()


_MARKER_RE = re.compile(r"__LG_JIT_KERNEL__\((\w+)\)")


def extract_jit_kernels(input_file):
    # Expand all macros, passing through any #include directives that can't be resolved
    expanded_code = _preprocess(input_file)

    lines = expanded_code.splitlines()
    kernels = {}

    i = 0
    while i < len(lines):
        m = _MARKER_RE.search(lines[i])
        if m:
            kernel_name = m.group(1)

            # Advance to the next non-empty line (the struct/function definition)
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1

            # Brace-count to find the end of the definition
            content = []
            brace_count = 0
            started = False

            while i < len(lines):
                line = lines[i]
                content.append(line)
                brace_count += line.count("{")
                brace_count -= line.count("}")

                if "{" in line:
                    started = True

                if started and brace_count == 0:
                    break
                i += 1

            # Join and format as a C string literal
            full_body = "\n".join(content)
            escaped_body = (
                full_body.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
            )
            kernels[kernel_name] = f'"{escaped_body}"'
        i += 1
    return kernels


def write_jit_header(input_file, output_file):
    """Generate a JIT header for input_file and write it to output_file."""
    kernels = extract_jit_kernels(input_file)
    with open(output_file, "w") as f:
        f.write("// This is a generated file containing all of the JIT strings\n")
        f.write("// for the kernels defined in " + input_file + "\n")
        f.write("#pragma once\n")
        for name, string_lit in kernels.items():
            f.write(f"static const char* {name}_JIT_STR = {string_lit};\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract LG_JIT_KERNEL-marked definitions from a C file as string literals."
    )
    parser.add_argument("--input", required=True, help="Input C source file")
    parser.add_argument(
        "--output", required=True, help="Output header file to generate"
    )
    args = parser.parse_args()
    write_jit_header(args.input, args.output)
