# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def InfoImporter(input_file):
    """
        Import a reference info file with the following format:
        #Chr name   NC   ID   Length
        chr1   NC_000001.9   Chr1#NC_000001.9#1#247249719#-1   247249719
        chr2   NC_000002.10   Chr2#NC_000002.10#1#242951149#-1   242951149

        returns a dict
        {'chr1': {'nc': 'NC_000001.9', 'length': }247249719 , ...}
    """

    reference_data = {}
    with open(input_file, 'r') as reference_info:
        header  = next(reference_info)
        header_map = {item[1].lower(): item[0] for item in enumerate(header.rstrip().lstrip("#").split("\t"))}
        for row in reference_info:
            columns = row.rstrip().split("\t")
            reference_data[columns[header_map['chr name']]] = {'nc': columns[header_map['nc']], 'length': int(columns[header_map['length']])}
    return reference_data
