import subprocess

# This is not particularly good code. Specifically, using shell=True is
# not suggested. That said, unless you're feeding insecure strings from
# the internet straight into this code, you're really fine. ACGU can't
# hurt you.

# You'll need to have Vienna installed and RNAfold on your path (indeed,
# RNAcofold for complex structures) to use this code.

def parse_score(seq2):
    ps = subprocess.check_output('echo "{}" | RNAfold'.format(seq2), shell=True).decode("utf-8")
    return float(ps.split('\n')[1].replace(' -', '-').replace('  0', '0').split()[1].replace('(', '').replace(')', ''))

def parse_structured_score(seq2, structure2):
    ps = subprocess.check_output('echo "{}\n{}" | RNAfold --enforceConstraint --constraint'.format(seq2, structure2.replace(">", "(").replace("<", ")")), shell=True).decode("utf-8")
    return float(ps.split('\n')[1].replace(' -', '-').replace('  0', '0').split()[1].replace('(', '').replace(')', ''))

def fold(seq2):
    ps = subprocess.check_output('echo "{}" | RNAfold'.format(seq2), shell=True).decode("utf-8")
    return ps.split('\n')[1].split()[0]

def parse_coscore(seq1, seq2):
    ps = subprocess.check_output('echo "{}&{}" | RNAcofold'.format(seq1, seq2), shell=True).decode("utf-8")

    secondline = ps.split('\n')[1].replace(' -', '-')
    secondline = secondline.split()[1].replace('(', '').replace(')', '')
    return float(secondline)
