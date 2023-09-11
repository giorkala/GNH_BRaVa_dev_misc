# GNH_BRaVa_dev_misc
Various scripts and pipelines developed for Genes &amp; Health as part of the BRaVa consortium.

In brief, the following work is undertaken:
1. Phasing of Genes & Health. This is performed using Shapeit5, following Lassen et al. (2023).
2. Variant annotation and encoding of additive/recessive genotypes. This is in accordance with [call_chets](https://github.com/frhl/call_chets). Annotation follows the [BRaVa guidelines](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit).
3. Testing for recessive effects. *Future work*

### Notes
* In `recessive_encoding/run_call_chets.sh`, the `cpp_dir/get_non_ref_sites` does almost the same as the proceeding line with bcftools, but the latter is safer to use and probably faster. 
* Occasionally, I've noted that some parts of the recessive encoding pipeline crash, giving an something crashes occasionally with throwing an "OSError: [Errno 12] Cannot allocate memory". Can't find why this happens, but re-running the corresponding task solves the problem.