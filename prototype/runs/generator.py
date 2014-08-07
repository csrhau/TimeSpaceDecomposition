#!/usr/bin/env python
from __future__ import print_function
import argparse
import errno
import os
import jinja2

def process_arguments():
    parser = argparse.ArgumentParser(description='QSUB Script Generator')
    parser.add_argument('--min_nodes', type=int, default=1)
    parser.add_argument('--max_nodes', type=int, required=True)
    parser.add_argument('--ppn', type=int, default=12)
    parser.add_argument('--cfg_tmpl', type=str, required=True)
    parser.add_argument('--pbs_tmpl', type=str, required=True)
    parser.add_argument('--config_dir', default='config')
    parser.add_argument('--pbs_dir', default='pbs')
    parser.add_argument('--name', type=str,
                        default='m_{{ mesh }}_n_{{ nodes }}')
    return parser.parse_args()

def generate_template(path):
    templateLoader = jinja2.FileSystemLoader(searchpath='./')
    j2env = jinja2.Environment(loader=templateLoader, trim_blocks=True)
    j2env.add_extension('jinja2.ext.loopcontrols')
    template = j2env.get_template(path)
    return template

def ensure_path_exists(path):
    # NOTE: redundant in python3, use  os.makedirs(path,exist_ok=True)
    try: 
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def main():
    args = process_arguments()
    cfg_tmpl = generate_template(args.cfg_tmpl)
    pbs_tmpl = generate_template(args.pbs_tmpl)
    name_tmpl = jinja2.Template(args.name)
    for nodes in xrange(args.min_nodes, args.max_nodes):
        for mesh in ['static', 'dynamic']:
            run_properties = {'mesh': mesh,
                              'nodes': nodes,
                              'ppn': args.ppn}
            name = name_tmpl.render(run_properties)
            run_properties['name'] = name
            ensure_path_exists(args.config_dir)
            ensure_path_exists(args.pbs_dir)
            cfg_filename = '{}/{}.cfg'.format(args.config_dir, name)
            pbs_filename = '{}/{}.pbs'.format(args.pbs_dir, name)
            with open(cfg_filename, 'w') as outfile:
                print(cfg_tmpl.render(run_properties), file=outfile)
            with open(pbs_filename, 'w') as outfile:
                print(pbs_tmpl.render(run_properties, config=cfg_filename),
                                      file=outfile)

if __name__ == '__main__':
  main()
