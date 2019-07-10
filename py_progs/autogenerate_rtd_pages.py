"""
This script goes through the created yaml documentation and converts them to
.rst files, then adds the new RST files to the git repository. It then calls
Sphinx to build the docs, in the documentation/html folder, and opens them in
the browser to test.

Arguments:
    Any
        Prints this documentation
"""
#!/usr/bin/env python
from typing import TextIO
import os
import sys
import webbrowser
from subprocess import call
import yaml


def write_header_by_level(output_file: TextIO, string: str, level: int = 0):
    """
    Writes a passed string as a header, with an appropriate type of underscore depending on
    nesting level (e.g. heading -> subheading)

    Arguments:
        output_file: The file to write to
        string: The text to write
        level: The level of heading, 0=heading, 1=subheading, 2=subsubheading, 3=paragraph heading
    """
    if level == 0:
        output_file.write("{}\n{}\n".format(string, ''.join('=' for i in range(len(string)))))
    elif level == 1:
        output_file.write("{}\n{}\n".format(string, ''.join('-' for i in range(len(string)))))
    elif level == 2:
        output_file.write("{}\n{}\n".format(string, ''.join('^' for i in range(len(string)))))
    else:
        output_file.write('**{}**\n{}\n'.format(string, ''.join('"' for i in range(len(string)+4))))


def write_str_indent(output_file: TextIO, string: str, indent: str = "  ", all: bool = False):
    """
    Writes a passed string to the output file. If it's split over multiple lines,
    indent the subsequent lines by the specified amount.

    Arguments:
        output_file: The file to write to
        string: The text to write
        indent: The indent to apply if the text splits over multiple lines
        all: Whether or not all lines should be indented
    """
    lines = string.splitlines()
    if all:
        output_file.write("{}{}\n".format(indent, lines[0]))
    else:
        output_file.write("{}\n".format(lines[0]))

    for line in lines[1:]:
        output_file.write("{}{}\n".format(indent, line))
    output_file.write("\n")
    return


def output_parameter(output_file: TextIO, parameter: dict, level: int = 0, no_separator: bool = False):
    """
    Output a parameter to file

    Arguments:
        output_file: The file to write to
        parameter: The parameter to write to file
        level: The level of the parameter, e.g. heading, subheading
    """
    # Suppress transition line at top level
    if not level and not no_separator:
        output_file.write("----------------------------------------\n\n")
    write_header_by_level(output_file, parameter['name'], level)
    output_file.write("{}\n".format(parameter['description']))

    if parameter.get('type'):
        output_file.write("**Type:** {}\n\n".format(parameter['type']))
    if parameter.get('unit'):
        output_file.write("**Unit:** {}\n\n".format(parameter['unit']))

    if parameter.get('values'):
        if isinstance(parameter['values'], dict):
            output_file.write("**Values:**\n\n")
            for key, value in parameter['values'].items():
                output_file.write("{}\n".format(key))
                if isinstance(value, str):
                    write_str_indent(output_file, value.strip(), indent="  ", all=True)
                elif isinstance(value, list):
                    write_str_indent(
                        output_file,
                        ', '.join([str(x) for x in value]), indent="  ", all=True
                    )
                else:
                    output_file.write("    "+str(value))

        elif isinstance(parameter['values'], list):
            # If this is a list of values, write each out as a bullet-point
            output_file.write("**Values:**\n")
            for value in parameter['values'].items():
                write_str_indent(output_file, "* {}".format(value), indent="  ")

        else:
            output_file.write("**Values:** {}\n".format(parameter['values']))
        output_file.write("\n")

    if parameter.get('parent'):
        if isinstance(parameter['parent'], dict):
            output_file.write("**Parent(s):**\n\n")
            for key, value in parameter['parent'].items():
                if isinstance(value, str):
                    write_str_indent(output_file, "* :ref:`{}`: {}".format(key, value), indent="  ")
                elif isinstance(value, list):
                    list_text = ', '.join(['``'+str(x)+'``' for x in value])
                    write_str_indent(output_file, "* :ref:`{}`: {}".format(key, list_text, indent="  "))
                else:
                    output_file.write("* :ref:`{}`: ``{}``\n\n".format(key, value))
            output_file.write("\n")

        elif isinstance(parameter['parent'], list):
            # If this is a list of parents, write each out as a bullet-point
            output_file.write("**Parent(s):**\n")
            for value in parameter['parent'].items():
                write_str_indent(output_file, "* {}".format(value), indent="  ")
        else:
            output_file.write("**Parent(s):** {}\n\n".format(parameter['parent']))

    output_file.write("**File:** {}\n\n\n".format(parameter['file']))

    # Go through all the children and output them
    if parameter.get('children'):
        for child in parameter['children'].values():
            output_parameter(output_file, child, level+1)
    return


def read_yaml(input_folder: str) -> dict:
    """
    Read all the .yaml files in to dicts. Categorise by type, e.g. 'reverb.mode' is type 'reverb'

    Arguments:
        input_folder: Path to the folder of input yaml files.

    Returns:
        dict: Dictionary of parameter dicts
    """
    dox_all = {}

    # Read in all parameters from the input files
    for input_file in [filename for filename in os.listdir(input_folder) if filename.endswith(".yaml")]:
        with open(os.path.join(input_folder, input_file), "r") as file_object:
            # Open each yaml file and load from it.
            try:
                parameter = yaml.load(file_object)
                parameter_name = parameter['name']

                # Get the type from the name (e.g. reverb.mode and set no children)
                parameter_type = parameter['name'].split('.')[0]
                parameter['root_type'] = parameter_type.strip()
                parameter['children'] = {}

                # Register this on the full list of all parameters
                dox_all[parameter_name] = parameter

            except yaml.YAMLError as error:
                print(error)

    return dox_all


def write_rst(output_folder: str, dox_structured: dict):
    """
    Write the structured documentation out to rst files

    Arguments:
        output_folder: The folder to write to
        dox_structured: The tree structured documentation
    """
    # Go through the parameter tree and write it out
    for parameter_type, type_parameters in dox_structured.items():
        with open(os.path.join(output_folder, "{}.autogen.rst".format(parameter_type)), 'w') as type_file:
            header_line = ''.join('=' for i in range(len(parameter_type)))
            type_file.write("{}\n{}\n{}\n\n".format(header_line, parameter_type, header_line))
            no_separator=True
            for parameter in type_parameters.values():
                output_parameter(type_file, parameter, level=0, no_separator=True)
                no_separator=False


def autogenerate_rtd_pages():
    """
    Write the RTD files to disk and add them to git, then run sphinx-build to generate the docs.
    """
    output_folder = os.path.join(os.environ["PYTHON"], "docs", "rst", "parameters")
    html_folder = os.path.join(os.environ["PYTHON"], "docs", "html")
    docs_folder = os.path.join(os.environ["PYTHON"], "docs")
    dox_all = {}

    try:
        # Try making output folder
        os.makedirs(output_folder)

    except OSError:
        # If we can't create the folder, clear out the old autogenerated files
        for filename in os.listdir(output_folder):
            if "autogen" in filename:
                os.remove(os.path.join(output_folder, filename))

        call(["git", "rm", "-f", os.path.join(output_folder, "*.autogen.rst")])

    dox_all = read_yaml(
        os.path.join(os.environ["PYTHON"], "docs", "parameters")
    )

    # This is a tree-structured dictionary e.g.
    # reverb: {
    #    children: {
    #        mode: {},
    #        visualisation: {
    #            children: {
    #                angular_bins: {}
    #            }
    #        }
    #    }
    # }
    dox_structured = {}

    # For each possible root type, set up a root dictionary for it
    for parameter_name, parameter in dox_all.items():
        dox_structured[parameter['root_type']] = {}

    # We want to build a tree with children etc. to respect parentage
    # Now, for the unstructured list of all parameters
    for parameter_name, parameter in dox_all.items():
        # For each parameter in the list of all parameters, check its parents
        parameter_type = parameter['root_type']

        if parameter.get('parent'):
            # If it has parents
            found_parent = False
            if isinstance(parameter['parent'], dict):
                for parent_name in parameter['parent'].keys():
                    # Check each parent.
                    if parent_name in dox_all.keys():
                        # If it's parent is another parameter (i.e. in .keys())
                        if parameter_type == dox_all[parent_name]['root_type']:
                            # If both the parameter and its parents are in the same class
                            # We assign this parameter to the children's list of its parent
                            dox_all[parent_name]["children"][parameter_name] = parameter
                            # And then 'continue' on, skipping other parents
                            found_parent = True
                            break

            if not found_parent:
                # If we didn't find the parameter's parent elsewhere,
                # then add it to the root.
                dox_structured[parameter_type][parameter_name] = parameter

        else:
            # If it has no parents, add it to the top layer
            dox_structured[parameter_type][parameter_name] = parameter

    write_rst(output_folder, dox_structured)

    call(["sphinx-build", "-b", "html", docs_folder, html_folder])
    call(["git", "add", os.path.join(output_folder, "*.autogen.rst")])
    webbrowser.open(os.path.join(html_folder, "index.html"))


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(__doc__)
    else:
        autogenerate_rtd_pages()
