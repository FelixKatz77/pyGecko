from ord_schema import message_helpers
import pytest

from pygecko.reaction import Reaction_Parser


def test_build_dataset_uses_current_ord_save_api(monkeypatch, tmp_path):
    dataset = object()
    output = tmp_path / 'dataset.pbtxt'
    writes = []
    monkeypatch.setattr(
        Reaction_Parser,
        'create_dataset_from_layout',
        classmethod(lambda cls, layout, yield_array: dataset))
    monkeypatch.setattr(
        message_helpers,
        'save_message',
        lambda message, path: writes.append((message, path)),
        raising=False)
    monkeypatch.setattr(
        message_helpers,
        'write_message',
        lambda message, path: pytest.fail('deprecated ORD writer was used'))

    result = Reaction_Parser.build_dataset(None, None, output)

    assert result is dataset
    assert writes == [(dataset, str(output))]
