import subprocess as sp


def run_librpa(args, cwd, out: str = "librpa.out", err: str = "librpa.err"):
    with open(cwd / out, "w") as stdo, open(cwd / err, "w") as stde:
        process = sp.Popen(args, cwd=cwd, stdout=stdo, stderr=stde)
        return_code = process.wait()
    return return_code
