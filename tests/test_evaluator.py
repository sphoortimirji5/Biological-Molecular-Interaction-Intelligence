"""Unit tests for ModelEvaluator — known predictions, deterministic assertions."""
import pytest
import numpy as np

from src.models.evaluator import ModelEvaluator


@pytest.fixture
def evaluator():
    return ModelEvaluator(k_values=[5, 10])


def test_perfect_predictions(evaluator):
    """Perfect classifier should score 1.0 everywhere."""
    y_true = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    y_prob = np.array([0.9, 0.8, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

    report = evaluator.evaluate(y_true, y_prob, split_name="test")

    assert report["threshold_free"]["auroc"] == 1.0
    assert report["threshold_free"]["auprc"] == 1.0
    assert report["threshold_based"]["accuracy"] == 1.0
    assert report["threshold_based"]["f1"] == 1.0


def test_random_predictions(evaluator):
    """Random predictions should score around 0.5 AUROC."""
    rng = np.random.RandomState(42)
    y_true = rng.randint(0, 2, 200)
    y_prob = rng.rand(200)

    report = evaluator.evaluate(y_true, y_prob, split_name="val")

    auroc = report["threshold_free"]["auroc"]
    assert 0.3 <= auroc <= 0.7, f"Random AUROC should be ~0.5, got {auroc}"


def test_ranking_metrics(evaluator):
    """Precision@k and recall@k with known ranking."""
    # Top 5 predictions: 4 positives, 1 negative
    y_true = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0])
    y_prob = np.array([0.95, 0.9, 0.85, 0.8, 0.75, 0.1, 0.1, 0.1, 0.1, 0.1])

    report = evaluator.evaluate(y_true, y_prob)
    ranking = report["ranking"]

    assert ranking["precision@5"] == 0.8    # 4/5 = 0.8
    assert ranking["recall@5"] == 1.0       # 4/4 = 1.0 (all positives found)
    assert ranking["precision@10"] == 0.4   # 4/10 = 0.4
    assert ranking["recall@10"] == 1.0      # 4/4 = 1.0


def test_report_structure(evaluator):
    """Report must contain all expected sections."""
    y_true = np.array([1, 0, 1, 0])
    y_prob = np.array([0.9, 0.1, 0.8, 0.2])

    report = evaluator.evaluate(y_true, y_prob, split_name="test")

    assert report["split"] == "test"
    assert report["n_samples"] == 4
    assert report["n_positive"] == 2
    assert report["n_negative"] == 2
    assert "threshold_free" in report
    assert "threshold_based" in report
    assert "ranking" in report
    assert "evaluated_at" in report


def test_single_class_returns_none(evaluator):
    """If only one class is present, AUROC/AUPRC should be None."""
    y_true = np.array([1, 1, 1, 1])
    y_prob = np.array([0.9, 0.8, 0.7, 0.6])

    report = evaluator.evaluate(y_true, y_prob)

    assert report["threshold_free"]["auroc"] is None
    assert report["threshold_free"]["auprc"] is None


def test_compare_reports(evaluator):
    """Compare produces side-by-side summary of multiple splits."""
    r1 = evaluator.evaluate(np.array([1, 0]), np.array([0.9, 0.1]), "train")
    r2 = evaluator.evaluate(np.array([1, 0]), np.array([0.8, 0.2]), "test")

    comparison = evaluator.compare([r1, r2])

    assert "train" in comparison["comparison"]
    assert "test" in comparison["comparison"]
    assert comparison["comparison"]["train"]["auroc"] is not None
    assert "compared_at" in comparison
