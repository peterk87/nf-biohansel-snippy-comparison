#!/usr/bin/env python
import click
import pandas as pd


def try_parse_int(x):
    try:
        return int(x)
    except ValueError:
        return x


def trace_entry(entry_string):
    out = {}
    for record_string in entry_string.strip().split('\n'):
        if record_string and '=' in record_string:
            key, value = record_string.split('=')
            key = key.strip()
            value = value.strip()
            value = try_parse_int(value)
            out[key] = value
    return out


def parse_trace_file(trace_file):
    with open(trace_file) as f:
        return [trace_entry(x) for x in f.read().strip().split('@@@') if x]

@click.command()
@click.option('-t', '--trace-file', required=True)
@click.option('-o', '--trace-out', required=True)
def main(trace_file, trace_out):
    traces = parse_trace_file(trace_file)
    df = pd.DataFrame(traces)
    df.to_csv(trace_out)

if __name__ == '__main__':
    main()
