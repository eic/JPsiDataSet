// Convert the first N events of a HepMC3 ROOT tree to HepMC3 ASCII.
// usage: tree2ascii IN.hepmc3.tree.root OUT.hepmc N
#include <HepMC3/GenEvent.h>
#include <HepMC3/ReaderRootTree.h>
#include <HepMC3/WriterAscii.h>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  if (argc < 4) { std::fprintf(stderr, "usage: %s IN OUT N\n", argv[0]); return 1; }
  HepMC3::ReaderRootTree in(argv[1]);
  HepMC3::WriterAscii out(argv[2], in.run_info());
  const long n = std::atol(argv[3]);
  long i = 0;
  for (; i < n; ++i) {
    HepMC3::GenEvent evt;
    if (!in.read_event(evt) || in.failed()) break;
    out.write_event(evt);
  }
  out.close();
  in.close();
  std::printf("wrote %ld events to %s\n", i, argv[2]);
}
