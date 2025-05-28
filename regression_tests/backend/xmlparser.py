from collections import OrderedDict
import xml.etree.ElementTree as ET

__all__ = ["XMLParser",]

TRUE_FALSE_NONE = {"true": True, "false": False, "none": None}


def _process_build_parameters(node_testcase):
    b = node_testcase.find('build')
    pars_tfn = ["require_libri", "require_greenx_api"]
    ret = {}
    for p in pars_tfn:
        try:
            value = b.get(p, None)
            value = TRUE_FALSE_NONE.get(value.lower(), value)
        except KeyError:
            value = None
        ret[p] = value
    return ret


def _process_run_parameters(node_testcase):
    r = node_testcase.find('run')
    ret = {}

    pars_exclud_intrange = ["ntasks_disable", "nthreads_disable"]
    for p in pars_exclud_intrange:
        try:
            value = r.get(p, None)
            if value.lower() in ["none", "false"]:
                value = tuple()
            else:
                value = tuple(int(x) for x in value.replace(',', ' ').split())
        except (KeyError, ValueError):
            value = None
        ret[p] = value
    return ret


def _process_validate_parameters(node_validate):
    nv = node_validate
    # name is manditory
    name = nv.get('name')
    file = nv.get('file', 'librpa.out')

    pars = ["regex", "occurence", "tablerows", "comparison", "binary"]
    ret = dict(name=name, file=file)
    for p in pars:
        ret[p] = nv.get(p, None)

    return ret


class XMLParser:

    def __init__(self, path_xml: str):
        tree = ET.parse(path_xml)
        root = tree.getroot()
        groups = OrderedDict()
        for group in root.findall('group'):
            gname = group.get('name')
            if gname not in groups:
                groups[gname] = []
            for tc in group.findall('testcase'):
                name = tc.get('name')
                directory = tc.get('directory')
                build = _process_build_parameters(tc)
                run = _process_run_parameters(tc)
                validates = []
                for v in tc.findall('validate'):
                    validates.append(_process_validate_parameters(v))
                groups[gname].append(dict(name=name, directory=directory, build=build, run=run, validates=validates))
        self.groups = groups

