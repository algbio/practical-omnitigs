from collections import OrderedDict
import itertools
import json
import traceback
import sys

# Safe formatting

class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'

def safe_format(str, **kwargs):
    return str.format_map(SafeDict(kwargs))

def safe_expand(str, **kwargs):
    items = []
    for key, values in kwargs.items():
        if type(values) is str or type(values) is not list:
            values = [values]
        items.append([(key, value) for value in values])

    for combination in itertools.product(*items):
        yield safe_format(str, **dict(combination))

def wildcard_format(str, wildcards):
    return str.format(**dict(wildcards.items()))

# Report file parsing

class Arguments(dict):
    @staticmethod
    def from_dict(arguments):
        assert type(arguments) is dict or type(arguments) is OrderedDict

        result = Arguments()
        for key, value in arguments.items():
            if type(value) is dict or type(value) is OrderedDict:
                result[key] = Arguments.from_dict(value)
            else:
                result[key] = value
        return result

    def update(self, other):
        #print("Calling update, other: {}".format(other))

        if other is None:
            return
        assert type(other) is Arguments

        SELECT_PREFIX = "select_"
        for key, value in other.items():
            if key.startswith(SELECT_PREFIX):
                key = key[len(SELECT_PREFIX):]
                select = True
            else:
                select = False

            if key not in self or self[key] is None:
                self[key] = value
            elif type(value) is str or type(value) is int or type(value) is float or type(value) is bool:
                # Selection logic: If only one value is given while we have a dict of values,
                # then select the dict entry with the respective key.
                if type(self[key]) is Arguments:
                    unselect = [k for k in self[key].keys() if k != value]
                    for k in unselect:
                        self[key].pop(k)

                    if len(self[key]) == 0:
                        sys.exit("Unselected all values: key: {}; unselect: {}".format(key, unselect))
                else:
                    assert type(self[key]) is str or type(self[key]) is int or type(self[key]) is float or type(self[key]) is bool, "type(self[key]) is not str, bool, int or float. key: {}, type(self[key]): {}".format(key, type(self[key]))
                    self[key] = value
            elif type(value) is dict or type(value) is OrderedDict or type(value) is Arguments or value is None:
                assert type(self[key]) is dict or type(self[key]) is OrderedDict or type(self[key]) is Arguments, f"wrong type(self[key]): {type(self[key])}, value type is: {type(value)}"
                
                # If in selection more, unselect all keys not part of the update.
                unselect = []
                if select:
                    if self[key] is None:
                        pass # TODO
                    unselect = [k for k in self[key].keys() if k not in value.keys()]
                    for k in unselect:
                        self[key].pop(k)

                self[key].update(value)

                if len(self[key]) == 0 and len(unselect) > 0:
                    sys.exit("Unselected all values: unselect: {}, selected_keys: {}".format(unselect, value.keys()))
            else:
                sys.exit("Cannot merge values of types {} and {}".format(type(self[key]), type(value)))

    def copy(self):
        result = Arguments()
        for key, value in self.items():
            if type(value) is str or type(value) is int or type(value) is float or type(value) is bool or value is None:
                result[key] = value
            elif type(value) is Arguments:
                result[key] = value.copy()
            else:
                sys.exit("Cannot copy value of type {}".format(type(value)))
        return result

    @staticmethod
    def argument_to_argument_string(key, value):
        if value is None or value == True:
            return str(key)
        elif value == False:
            return ""
        else:
            return str(key) + " " + str(value)

    def to_argument_string(self):
        result = ""
        once = True
        for key, value in self.items():
            argument = Arguments.argument_to_argument_string(key, value)
            if argument != "":
                if once:
                    once = False
                else:
                    result += " "
                result += argument
        return result

    def subarguments_name(self, name):
        if name not in self:
            return None

        result = self[name]
        if type(result) is Arguments:
            result = list(result.keys())
            if len(result) != 1:
                return None
            result = result[0]
        return result

    def subarguments_arguments(self, name):
        if name not in self:
            return None

        result = self[name]
        if type(result) is not Arguments:
            return Arguments()

        result = list(result.values())
        if len(result) != 1 or type(result[0]) is not Arguments:
            return Arguments()

        return result[0]

    def read_simulator_name(self):
        return self.subarguments_name("read_simulator")

    def read_simulator_arguments(self):
        return self.subarguments_arguments("read_simulator")

    def assembler_name(self):
        return self.subarguments_name("assembler")

    def assembler_arguments(self):
        return self.subarguments_arguments("assembler")

    def postprocessor_name(self):
        return self.subarguments_name("postprocessor")

    def postprocessor_arguments(self):
        return self.subarguments_arguments("postprocessor")

    def genome(self):
        if "genome" in self:
            return self["genome"]
        else:
            return None

    def homopolymer_compression(self):
        if "homopolymer_compression" in self:
            return self["homopolymer_compression"] == True # Compare with True to ensure return type is boolean
        else:
            return False

    def __str__(self):
        if self is None:
            return "{}"

        string = json.dumps(self, separators = (',', ':'))
        n = 200
        result = "/_/".join([string[i:i+n] for i in range(0, len(string), n)])
        return result

    def from_str(string):
        if string is None:
            return None
        elif string == "None":
            return None

        string = string.replace("/_/", "")
        return Arguments.from_dict(json.loads(string))

    def shortstr(self):
        string = ""
        shortened = self.copy()
        shortened.pop("genome", None)
        cli_arguments = shortened.get("cli_arguments", None)
        if type(cli_arguments) is Arguments:
            for k, v in cli_arguments.items():
                shortened[k] = v
            shortened.pop("cli_arguments")

        for key, value in shortened.items():
            if len(string) > 0:
                string += ","
            string += key.replace("-", "")[0] + ":"

            if type(value) is Arguments:
                string += "(" + value.shortstr() + ")"
            elif type(value) is int or type(value) is float:
                string += str(value)
            elif type(value) is bool:
                assert value == True
                string = string[:-1]
            else:
                string += str(value)[0:2]
        return string

    def retain_raw_assembly_arguments(self):
        if "postprocessor" in self:
            self.pop("postprocessor")

class Algorithm:
    def __init__(self, assembler, arguments):
        self.assembler = assembler
        self.arguments = arguments

class Column:
    def __init__(self, root_arguments, additional_arguments):
        #print("Creating column\nroot: {}\nadditional: {}".format(root_arguments, additional_arguments))
        assert type(root_arguments) is Arguments
        assert type(additional_arguments) is Arguments
        self.arguments = root_arguments.copy()
        self.arguments.update(additional_arguments.copy())
        #print("Result arguments: {}".format(self.arguments))
        self.shortname = self.arguments.pop("shortname")
        self.assembler = self.arguments.assembler_name()
        self.assembler_arguments = self.arguments.assembler_arguments()
        self.genome = self.arguments["genome"]

    def __str__(self):
        return str(self.arguments)

    def from_str(string):
        return Column(Arguments.from_str(string), Arguments())

class ReportFile:
    def __init__(self, arguments, columns):
        self.arguments = arguments
        self.columns = columns
        self.name = str(arguments)
        self.shortname = arguments.shortstr()

    def __str__(self):
        return "ReportFile(name: {}, arguments: {}, columns: {})".format(self.name, self.arguments, self.columns)

class ArgumentMatrix:
    def __init__(self, argument_matrix):
        self.argument_matrix = argument_matrix
        self.len = ArgumentMatrix._compute_length(argument_matrix)

    def __len__(self):
        return self.len

    def _compute_length(argument_matrix):
        if argument_matrix is None:
            return 1

        if type(argument_matrix) is ArgumentMatrix:
            argument_matrix = argument_matrix.argument_matrix

        if type(argument_matrix) is dict or type(argument_matrix) is OrderedDict:
            length = 1
            for key, value in argument_matrix.items():
                length *= ArgumentMatrix._compute_length(value)
            return length

        if type(argument_matrix) is list:
            length = 0
            for value in argument_matrix:
                length += ArgumentMatrix._compute_length(value)
            return length

        if type(argument_matrix) is str or type(argument_matrix) is int or type(argument_matrix) is float or type(argument_matrix) is bool:
            return 1

        sys.exit("Illegal type of argument matrix: {}.\nAllowed are dict, OrderedDict, list, str, bool, int and float".format(type(argument_matrix)))

    def __iter__(self):
        return ArgumentMatrix._iter(self.argument_matrix)

    def _iter(argument_matrix):
        if argument_matrix is None:
            return None

        if type(argument_matrix) is ArgumentMatrix:
            argument_matrix = argument_matrix.argument_matrix

        if type(argument_matrix) is str or type(argument_matrix) is int or type(argument_matrix) is float or type(argument_matrix) is bool:
            #print("Found str, int or float: {}".format(argument_matrix))
            yield argument_matrix

        elif type(argument_matrix) is dict or type(argument_matrix) is OrderedDict or type(argument_matrix) is Arguments:
            #print("Found dict: {}".format(argument_matrix))
            keys = [None] * len(argument_matrix)
            values = [None] * len(argument_matrix)
            for index, (key, value) in enumerate(argument_matrix.items()):
                keys[index] = key
                if value is None:
                    values[index] = [None]
                else:
                    values[index] = list(ArgumentMatrix._iter(value))

            #print("keys: {}".format(keys))
            #print("values: {}".format(values))

            for value_combination in itertools.product(*values):
                result = Arguments()
                for key, value in zip(keys, value_combination):
                    if value is not None:
                        result[key] = value

                if len(result) == 1:
                    if list(result.values())[0] is None:
                        yield list(result.keys())[0]
                        continue

                yield result

        elif type(argument_matrix) is list:
            #print("Found list: {}".format(argument_matrix))
            for value in argument_matrix:
                if value is None:
                    yield None
                else:
                    for list_item in ArgumentMatrix._iter(value):
                        yield list_item

        else:
            sys.exit("Illegal type of argument matrix: {}.\nAllowed are dict, OrderedDict, list, str, bool, int and float".format(type(argument_matrix)))