"""
Run with:
    pip install -r requirements.txt
    pytest test_external_usage.py -v
"""
from skipprobabilities import run_pipeline, DerivationPipeline

def test_import_works_from_outside_the_repo():
    assert run_pipeline is not None
    assert DerivationPipeline is not None

def test_pipeline_runs_from_outside_the_repo():
    derivation = run_pipeline(output_dir="./external_test_out")
    assert isinstance(derivation, DerivationPipeline)
    assert len(derivation.skip_probs) > 0