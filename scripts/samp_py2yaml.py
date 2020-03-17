#!/usr/bin/env python

"""
Helper script for converting sample+samplegroup loader scripts into YAML files specifying which samples+samplegroups to load

Args:
    --inpy: python file that loaded all the samples/sample groups
    --outpfx: prefix of the output yml files (i.e. ./res/r17 -> ./res/r17_samples.yml and ./res/r17_samplegroups.yml)

Author: Ryan Tse
"""

import re
import argparse

import yaml


parser = argparse.ArgumentParser()
parser.add_argument('--inpy', help='the python file that loaded all the samples/sample groups')
parser.add_argument('--outpfx', help='prefix of the output yml files (i.e. ./res/r17 -> ./res/r17_samples.yml and ./res/r17_samplegroups.yml)')
args = parser.parse_args()


def split_exclude(instr, spstr, exstrs):
    """
    Split an input string (`instr`) into component parts, breaking along instances of `spstr` except for cases
    where `spstr` occurs between delimiters (i.e. parenthesis or quote marks) specified as a list in exstrs
    
    Ex.
        split_exclude("hello, (there, bob) how, are you, $doing, I wonder$", ",", ["()", "$$"])
        = ["hello", " (there, bob) how", " are you", "$doing, I wonder$"]
    """
    out = []
    rep = ''
    within = False
    delim = ""
    for char in instr:
        if within:
            if char == delim[1]:
                within = False
            rep += char
            continue
            
        for i in range(len(exstrs)):
            if char == exstrs[i][0]:
                delim = exstrs[i]
                within = True
                continue
        
        if char == spstr:
            out.append(rep)
            rep = ''
            continue
        
        rep += char
    out.append(rep)
    
    return out


def main():
    with open(args.inpy, 'r') as fi:
        fs = fi.read()

        # Strip comments
        b = 0
        fsr = ''
        for comment in re.finditer("#(?=([^']*'[^']*')*[^']*$)", fs[b:]):
            ci = comment.start(0)
            fsr += fs[b:ci]
            b = fs.find('\n', ci)
        fsr += fs[b:]
        fs = fsr

        # Ignore up to the first added sample
        beg = fs.find("samples.AddSample(") + len("samples.AddSample(")

        # Split all the AddSample
        entries = fs[beg:fs.find('samples.AddSampleGroup')].split("samples.AddSample(")    
        samples = []
        for entry in entries:
            entry = entry.strip()
            entry = entry[:-1]    # trailing paren

            # Clean up whitespace       
            rep = split_exclude(entry, ",", ['""', "''", '[]'])
            rep = [''.join(rep_s.split()) for rep_s in rep]

            # Assemble dict
            sample = {'name': eval(rep[0])}
            for kv in rep[1:]:
                kvl = kv.split('=')
                kvl = [kvls.strip() for kvls in kvl]

                if kvl[1][0] == '"' or kvl[1][0] == "'":
                        kvl[1] = kvl[1][1:-1]

                kvl_copy = '%s' % kvl[1]    # make a (deep) copy
                try:
                    kvl[1] = eval(kvl[1])
                except:
                    kvl[1] = kvl_copy

                sample[kvl[0]] = kvl[1]

            samples.append(sample)

        # Split all the AddSampleGroup
        beg = fs.find("samples.AddSampleGroup(") + len("samples.AddSampleGroup(")
        entries = fs[beg:fs.find("def print_examples() :")].split('samples.AddSampleGroup(')    
        sample_groups = []
        for entry in entries:        
            entry = entry.strip()
            entry = entry[:-1]    # trailing paren

            # Clean up whitespace
            rep = split_exclude(entry, ",", ['""', "''", '[]'])
            rep = [''.join(rep_s.split()) for rep_s in rep]

            # Assemble dict
            sample_group = {'name': eval(rep[0])}
            for kv in rep[1:]:
                if kv != '':
                    kvl = kv.split('=')
                    kvl = [kvls.strip() for kvls in kvl]

                    if kvl[1][0] == '"' or kvl[1][0] == "'":
                        kvl[1] = kvl[1][1:-1]

                    kvl_copy = '%s' % kvl[1]    # make a (deep) copy
                    try:
                        kvl[1] = eval(kvl[1])
                    except:
                        kvl[1] = kvl_copy

                    sample_group[kvl[0]] = kvl[1]

            sample_groups.append(sample_group)

        # Save to file
        with open(args.outpfx + '_samples.yml', 'w') as of:
            yaml.dump(samples, of, default_flow_style=False)
        with open(args.outpfx + '_samplegroups.yml', 'w') as of:
            yaml.dump(sample_groups, of, default_flow_style=False)

if __name__ == "__main__":
    main()