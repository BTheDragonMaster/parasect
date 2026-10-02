# -*- coding: utf-8 -*-

"""Unit tests for parasect.core.helpers: prepare_folders, prepare_model, and
download_and_unpack_or_fetch (used by cli_paras.py to set up the temp/model
directories and fetch/retrain the PARAS model before running predictions).

Out of scope: prepare_substrates() -- it exists in this module but is only used
by the PARASECT custom-substrate path, not by `paras`.
"""

import logging
import os
import tarfile
import unittest
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, patch

from parasect.core.helpers import download_and_unpack_or_fetch, prepare_folders, prepare_model
from parasect.core.models import Model, ModelType

logger = logging.getLogger("test_helpers")


class TestPrepareFolders(unittest.TestCase):
    def test_creates_out_dir_and_default_temp_and_model_subdirs(self):
        with TemporaryDirectory() as tmp_dir:
            out_dir = os.path.join(tmp_dir, "output")
            temp_dir, model_dir = prepare_folders(out_dir, None, None)

            self.assertEqual(temp_dir, os.path.join(out_dir, "temp"))
            self.assertEqual(model_dir, os.path.join(out_dir, "model"))
            self.assertTrue(os.path.isdir(out_dir))
            self.assertTrue(os.path.isdir(temp_dir))
            self.assertTrue(os.path.isdir(model_dir))

    def test_uses_given_temp_and_model_dirs_when_provided(self):
        with TemporaryDirectory() as tmp_dir:
            out_dir = os.path.join(tmp_dir, "output")
            given_temp = os.path.join(tmp_dir, "my_temp")
            given_model = os.path.join(tmp_dir, "my_model")

            temp_dir, model_dir = prepare_folders(out_dir, given_temp, given_model)

            self.assertEqual(temp_dir, given_temp)
            self.assertEqual(model_dir, given_model)
            self.assertTrue(os.path.isdir(given_temp))
            self.assertTrue(os.path.isdir(given_model))

    def test_is_idempotent_when_directories_already_exist(self):
        with TemporaryDirectory() as tmp_dir:
            out_dir = os.path.join(tmp_dir, "output")
            prepare_folders(out_dir, None, None)
            # calling again with the same, now-existing directories must not raise
            temp_dir, model_dir = prepare_folders(out_dir, None, None)
            self.assertTrue(os.path.isdir(temp_dir))
            self.assertTrue(os.path.isdir(model_dir))

    def test_raises_filenotfounderror_if_given_temp_dirs_parent_is_missing(self):
        # os.mkdir (not os.makedirs) is used internally, so a multi-level path
        # whose parent doesn't exist yet will fail rather than being created.
        # Documenting this real limitation rather than changing it.
        with TemporaryDirectory() as tmp_dir:
            out_dir = os.path.join(tmp_dir, "output")
            nested_temp = os.path.join(tmp_dir, "does_not_exist_yet", "temp")
            with self.assertRaises(FileNotFoundError):
                prepare_folders(out_dir, nested_temp, None)


class TestDownloadAndUnpackOrFetch(unittest.TestCase):
    @patch("parasect.core.helpers.urlopen")
    def test_returns_existing_plain_file_without_downloading(self, mock_urlopen):
        mock_urlopen.side_effect = AssertionError("should not attempt to download an already-present file")
        with TemporaryDirectory() as tmp_dir:
            existing_file = os.path.join(tmp_dir, "model.paras.gz")
            with open(existing_file, "wb") as fo:
                fo.write(b"already here")

            result = download_and_unpack_or_fetch("https://example.com/model.paras.gz", tmp_dir, logger)

            self.assertEqual(result, existing_file)
            mock_urlopen.assert_not_called()

    @patch("parasect.core.helpers.urlopen")
    def test_returns_existing_extracted_archive_without_redownloading(self, mock_urlopen):
        mock_urlopen.side_effect = AssertionError("should not attempt to download an already-extracted archive")
        with TemporaryDirectory() as tmp_dir:
            # build a real tar.gz with a single root folder, plus that folder
            # already "extracted" alongside it, mirroring a prior successful run.
            archive_path = os.path.join(tmp_dir, "bundle.tar.gz")
            extracted_dir = os.path.join(tmp_dir, "bundle_root")
            os.mkdir(extracted_dir)
            with open(os.path.join(extracted_dir, "inner.txt"), "w") as fo:
                fo.write("hello")

            with tarfile.open(archive_path, "w:gz") as tar:
                tar.add(extracted_dir, arcname="bundle_root")

            result = download_and_unpack_or_fetch("https://example.com/bundle.tar.gz", tmp_dir, logger)

            self.assertEqual(result, extracted_dir)
            mock_urlopen.assert_not_called()

    def test_downloads_and_extracts_new_archive_with_single_root(self):
        with TemporaryDirectory() as source_dir, TemporaryDirectory() as dest_dir:
            # Build a real tar.gz to serve as the "remote" content.
            payload_dir = os.path.join(source_dir, "bundle_root")
            os.mkdir(payload_dir)
            with open(os.path.join(payload_dir, "inner.txt"), "w") as fo:
                fo.write("hello from archive")
            archive_path = os.path.join(source_dir, "bundle.tar.gz")
            with tarfile.open(archive_path, "w:gz") as tar:
                tar.add(payload_dir, arcname="bundle_root")

            with open(archive_path, "rb") as fo:
                archive_bytes = fo.read()

            mock_response = MagicMock()
            mock_response.__enter__.return_value = mock_response
            mock_response.getheader.return_value = str(len(archive_bytes))
            # first read() returns the full payload, second returns b"" to end the loop
            mock_response.read.side_effect = [archive_bytes, b""]

            with patch("parasect.core.helpers.urlopen", return_value=mock_response) as mock_urlopen:
                result = download_and_unpack_or_fetch("https://example.com/bundle.tar.gz", dest_dir, logger)

            mock_urlopen.assert_called_once()
            self.assertEqual(result, os.path.join(dest_dir, "bundle_root"))
            with open(os.path.join(result, "inner.txt")) as fo:
                self.assertEqual(fo.read(), "hello from archive")

    def test_downloads_new_plain_file_when_not_already_present(self):
        with TemporaryDirectory() as dest_dir:
            file_bytes = b"model contents"
            mock_response = MagicMock()
            mock_response.__enter__.return_value = mock_response
            mock_response.getheader.return_value = str(len(file_bytes))
            mock_response.read.side_effect = [file_bytes, b""]

            with patch("parasect.core.helpers.urlopen", return_value=mock_response) as mock_urlopen:
                result = download_and_unpack_or_fetch("https://example.com/model.paras.gz", dest_dir, logger)

            mock_urlopen.assert_called_once()
            self.assertEqual(result, os.path.join(dest_dir, "model.paras.gz"))
            with open(result, "rb") as fo:
                self.assertEqual(fo.read(), file_bytes)


class TestPrepareModel(unittest.TestCase):
    def test_raises_for_composite_model_type_before_touching_network(self):
        # model_needs_retraining() rejects composite ModelTypes up front.
        with TemporaryDirectory() as model_dir:
            with self.assertRaises(ValueError):
                prepare_model(ModelType.ANTISMASH_MODELS, model_dir)

    @patch("parasect.core.helpers.download_and_unpack_or_fetch")
    @patch("parasect.core.helpers.model_needs_retraining", return_value=False)
    def test_copies_metadata_file_into_model_dir_when_absent(self, mock_needs_retraining, mock_download):
        mock_download.return_value = "/fake/model/path"
        with TemporaryDirectory() as model_dir:
            self.assertFalse(os.path.exists(os.path.join(model_dir, "model_metadata.txt")))
            prepare_model(ModelType.PARAS, model_dir)
            self.assertTrue(os.path.exists(os.path.join(model_dir, "model_metadata.txt")))

    @patch("parasect.core.helpers.download_and_unpack_or_fetch")
    @patch("parasect.core.helpers.model_needs_retraining", return_value=False)
    def test_downloads_paras_model_from_expected_zenodo_url(self, mock_needs_retraining, mock_download):
        mock_download.return_value = "/fake/model/path"
        with TemporaryDirectory() as model_dir:
            result = prepare_model(ModelType.PARAS, model_dir)

            self.assertEqual(result, "/fake/model/path")
            called_url = mock_download.call_args[0][0]
            self.assertIn("model.paras.gz", called_url)
            self.assertNotIn("all_substrates", called_url)

    @patch("parasect.core.helpers.download_and_unpack_or_fetch")
    @patch("parasect.core.helpers.model_needs_retraining", return_value=False)
    def test_downloads_paras_all_substrates_model_from_expected_zenodo_url(self, mock_needs_retraining, mock_download):
        mock_download.return_value = "/fake/model/path"
        with TemporaryDirectory() as model_dir:
            prepare_model(ModelType.PARAS_ALL_SUBSTRATES, model_dir)

            called_url = mock_download.call_args[0][0]
            self.assertIn("all_substrates_model.paras.gz", called_url)

    @patch("parasect.core.helpers.download_and_unpack_or_fetch")
    @patch("parasect.core.helpers.model_needs_retraining", return_value=False)
    def test_raises_for_unrecognized_single_flag_model_type(self, mock_needs_retraining, mock_download):
        # A single-bit ModelType that isn't one of the four handled branches
        # (e.g. a future model type added to the enum but not to this function)
        # should be rejected rather than silently downloading nothing.
        unrecognized_model_type = ModelType(16)
        with TemporaryDirectory() as model_dir:
            with self.assertRaises(ValueError):
                prepare_model(unrecognized_model_type, model_dir)
        mock_download.assert_not_called()

    @patch("parasect.core.helpers.update_metadata_file")
    @patch("parasect.core.helpers.retrain_model")
    @patch("parasect.core.helpers.model_needs_retraining", return_value=True)
    def test_retrains_and_updates_metadata_when_retraining_is_needed(
        self, mock_needs_retraining, mock_retrain_model, mock_update_metadata
    ):
        fake_model = MagicMock(spec=Model)
        fake_model.file_name = "model.paras.gz"
        mock_retrain_model.return_value = fake_model

        with TemporaryDirectory() as model_dir:
            result = prepare_model(ModelType.PARAS, model_dir)

            mock_retrain_model.assert_called_once_with(ModelType.PARAS)
            fake_model.save.assert_called_once_with(model_dir)
            mock_update_metadata.assert_called_once()
            self.assertEqual(result, os.path.join(model_dir, "model.paras.gz"))


if __name__ == "__main__":
    unittest.main()
