#!/usr/bin/env python2
from __future__ import print_function
import os
import sys
import trace
import atexit
import tempfile
import subprocess
from shutil import copy, rmtree
from collections import defaultdict

# only needed to locate CMSSW
import six
import FWCore
import FWCore.ParameterSet.Types
import FWCore.Modules.EmptySource_cfi

OUTFILE_TREE = "calltree"
OUTFILE_FILES = "callfiles"
WRAP_SCRIPTS = ["cmsDriver.py" ]
IGNORE_DIRS = [
  os.path.dirname(os.__file__),
  os.path.dirname(six.__file__),
  FWCore.ParameterSet.Types.__file__,
]
STRIPPATHS = [
  "/".join(os.path.dirname(FWCore.__file__).split("/")[:-1]) + "/",
  "/".join(FWCore.Modules.EmptySource_cfi.__file__.split("/")[:-3]) + "/",
]


def setupenv():
  bindir = tempfile.mkdtemp()
  print("+Setting up in ", bindir)
  for s in WRAP_SCRIPTS:
    os.symlink(__file__, bindir + "/" + s)
  os.symlink(__file__, bindir + "/cmsRun")
  os.environ["PATH"] = bindir + ":" + os.environ["PATH"]
  os.environ["CMSSWCALLTREE"] = bindir + "/" + OUTFILE_TREE
  os.environ["CMSSWCALLFILES"] = bindir + "/" + OUTFILE_FILES
  with open(os.environ["CMSSWCALLTREE"], "w") as f:
    pass
  with open(os.environ["CMSSWCALLFILES"], "w") as f:
    pass
  return bindir

def cleanupenv(tmpdir):
  #with open(os.environ["CMSSWCALLTREE"], "a") as f:
  #  print("}", file=f)
  print("+Cleaning up ", tmpdir)
  copy(os.environ["CMSSWCALLTREE"], ".")
  copy(os.environ["CMSSWCALLFILES"], ".")
  rmtree(tmpdir)


def trace_command(argv):
  tmpdir = None
  if not "CMSSWCALLTREE" in os.environ:
    tmpdir = setupenv()

  subprocess.call(argv)

  if tmpdir:
    cleanupenv(tmpdir)

def funcid(frame):
  code = frame.f_code
  filename = code.co_filename
  funcname = code.co_name
  # compared to the `trace` module, we don't attempt to find class names here 
  return (filename, funcname)

def trace_python(prog_argv, path):
  files = set()
  callgraph = defaultdict(lambda: set())

  def nop_trace(frame, why, arg):
    pass

  def tracefunc(frame, why, arg):
    if why == 'call':
      code = frame.f_code
      # compared to the `trace` module, we don't attempt to find class names here 
      filename = code.co_filename

      for d in IGNORE_DIRS:
        if filename.startswith(d):
          sys.settrace(nop_trace)
          return wait_for_return

      funcname = code.co_name
      code = frame.f_back.f_code
      p_filename = code.co_filename
      p_funcname = code.co_name

      files.add(filename)
      callgraph[(filename, funcname)].add((p_filename, p_funcname))
    return None

  def wait_for_return(frame, why, arg):
    if why == 'return':
      sys.settrace(tracefunc)
    return wait_for_return

  sys.argv = prog_argv
  progname = prog_argv[0]

  # Search $PATH. There seems to be no pre-made function for this.
  for entry in path:
    file_path = os.path.join(entry, progname)
    if os.path.isfile(file_path):
      break
  if not os.path.isfile(file_path):
    print("+Cannot find program (%s) in modified $PATH (%s)." % (progname, path))
    sys.exit(1)
  print("+Found %s as %s in %s." % (progname, file_path, path))

  def writeoutput():
    print("+Done running %s, writing output..." % progname)

    def formatfile(filename):
      filename = os.path.abspath(filename)
      for pfx in STRIPPATHS:
        if filename.startswith(pfx):
          filename = filename[len(pfx):]
      return filename
      
    def format(func):
      filename, funcname = func
      return "%s::%s" % (formatfile(filename), funcname)

    def callpath(func):
      # climb up in the call graph until we find a node without callers (this is
      # the entry point, the traced call itself). There may be cycles, but any
      # node is reachable from the entry point, so no backtracking required.
      path = []
      seen = set()
      parents = {func}
      while parents:
        if len(parents) == 1:
          func = next(iter(parents))
          seen.add(func)
          path.append(format(func))
        if len(parents) > 1:
          for func in parents:
            if not func in seen:
              break
          seen.add(func)
          path.append(format(func) + "+")
        parents = callgraph[func]
      return path[:-1]

    with open(os.environ["CMSSWCALLFILES"], "a") as outfile:
        for f in files:
          print("%s: %s" % (progname, formatfile(f)), file=outfile)
    with open(os.environ["CMSSWCALLTREE"], "a") as outfile:
        for func in callgraph.keys():
          print("%s: %s 1" % (progname, ";".join(reversed(callpath(func)))), file=outfile)

  try:
    with open(file_path) as fp:
      code = compile(fp.read(), progname, 'exec')
      # try to emulate __main__ namespace as much as possible
      globals = {
      '__file__': progname,
      '__name__': '__main__',
      '__package__': None,
      '__cached__': None,
      }

      # would be too easy if this covered all the cases...
      atexit.register(writeoutput)
      # cmsDriver calls cmsRun via exec (execvpe specifically), so we also need
      # to hook that...
      old_execvpe = os.execvpe
      def exec_hook(*args):
        writeoutput()
        old_execvpe(*args)
      os.execvpe = exec_hook

      # now turn on the traceing
      sys.settrace(tracefunc)
      try:
        exec code in globals, globals
      finally:
        sys.settrace(None)

  except OSError as err:
    print("+Cannot run file %r because: %s" % (sys.argv[0], err))
    sys.exit(1)
  except SystemExit:
    pass
  # this is not necessarily reached at all. 
  sys.exit(0)

def help():
  print("Usage: %s <some cmssw commandline>" % (sys.argv[0])))
  print("  The given programs will be executed, instrumenting calls to %s and cmsRun." % (", ".join(WRAP_SCRIPTS)))
  print("  cmsRun will not actually run cmssw, but all the Python code will be executed and instrumentd. The results are written to the file `%s` in the same directory." % OUTFILE)
  print("Examples:")
  print("  %s --files runTheMatrix.py -l 1000 --ibeos" % sys.argv[0])
  print(  "%s --callgraph cmsRun rpc_dqm_sourceclient-live_cfg.py" % sys.argv[0])

def main():
  print("+Running cmsswfiletrace...")
  for s in WRAP_SCRIPTS:
    if sys.argv[0].endswith(s):
      print("+Wrapping %s..." % s)
      tmppath = os.path.dirname(sys.argv[0])
      path = filter(
        lambda s: not s.startswith(tmppath),
        os.environ["PATH"].split(":")
      )
      trace_python([s] + sys.argv[1:], path)
      return
  if sys.argv[0].endswith('cmsRun'):
      print("+Wrapping cmsRun...")
      trace_python(sys.argv[1:], ["."])
      return
  if len(sys.argv) <= 1:
    help()
    return
  # else
  print("+Running command with tracing %s..." % sys.argv[1:])
  trace_command(sys.argv[1:])


if __name__ == '__main__':
  main()
  
