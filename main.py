import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from external_tests.run_example import run_pipeline

derivation = run_pipeline(output_dir="./results/my_output")
print(derivation.print_blinded())
derivation.stats()


print(derivation.skip_probs)
print(derivation.skips)