"""Unit tests for FeatureStore — mock storage, no Docker needed."""
import io
import json
import pytest
import pandas as pd
from unittest.mock import MagicMock, ANY

from src.features.store import FeatureStore


@pytest.fixture
def mock_storage():
    return MagicMock()


@pytest.fixture
def feature_store(mock_storage):
    return FeatureStore(storage=mock_storage, bucket="test-features")


def test_save_features_and_manifest(feature_store, mock_storage):
    df = pd.DataFrame({"id": [1, 2], "value": [0.1, 0.2]})
    manifest = {"param": "value"}

    uri = feature_store.save("compound", "v1", df, manifest)

    assert uri == "s3://test-features/compound/v1/features.parquet"
    assert mock_storage.put_object.call_count == 2

    mock_storage.put_object.assert_any_call(
        bucket="test-features",
        key="compound/v1/features.parquet",
        body=ANY,
        content_type="application/vnd.apache.parquet",
    )
    mock_storage.put_object.assert_any_call(
        bucket="test-features",
        key="compound/v1/manifest.json",
        body=ANY,
        content_type="application/json",
    )


def test_load_features(feature_store, mock_storage):
    df_orig = pd.DataFrame({"id": [1, 2], "value": [0.1, 0.2]})
    buffer = io.BytesIO()
    df_orig.to_parquet(buffer)
    buffer.seek(0)

    mock_storage.get_object.return_value = buffer
    df_loaded = feature_store.load("compound", "v1")

    pd.testing.assert_frame_equal(df_orig, df_loaded)
    mock_storage.get_object.assert_called_with(
        bucket="test-features", key="compound/v1/features.parquet"
    )


def test_load_manifest(feature_store, mock_storage):
    manifest_orig = {"param": "value"}
    mock_storage.get_object.return_value = io.BytesIO(
        json.dumps(manifest_orig).encode("utf-8")
    )

    manifest_loaded = feature_store.load_manifest("compound", "v1")

    assert manifest_loaded == manifest_orig
    mock_storage.get_object.assert_called_with(
        bucket="test-features", key="compound/v1/manifest.json"
    )
