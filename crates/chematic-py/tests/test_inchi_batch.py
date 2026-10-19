"""Retain row positions, exact identities and typed failures for InChI batches."""

from dataclasses import FrozenInstanceError
import numpy as np
import pytest
import chematic as ch
import chematic._inchi_batch as implementation

WATER = "InChI=1S/H2O/h1H2"
ETHANOL = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"


def test_mixed_rows_keep_positions_sources_and_duplicates():
    sources = [WATER, "not-an-inchi", ETHANOL, "", WATER]
    rows = ch.inchi_keys_checked(iter(sources))
    assert isinstance(rows, tuple)
    assert [row.index for row in rows] == list(range(5))
    assert [row.source_text for row in rows] == sources
    assert [row.status for row in rows] == ["ok", "invalid", "ok", "invalid", "ok"]
    assert [row.key for row in rows] == [
        "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        None,
        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        None,
        "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
    ]
    assert rows[0].code is None and rows[0].message is None
    assert rows[1].code == "inchi_prefix"
    with pytest.raises(FrozenInstanceError):
        rows[0].key = "changed"


@pytest.mark.parametrize(
    "text,status,code",
    [
        ("InChI=1/CH4/h1H4", "unsupported", "standard_inchi_required"),
        (WATER + "\n", "invalid", "invalid_characters"),
        (" " + WATER, "invalid", "inchi_prefix"),
        ("\ud800", "invalid", "invalid_utf8"),
        ("InChI=1S/", "invalid", "invalid_formula"),
        ("InChI=1S/[CH4]", "invalid", "invalid_formula"),
        ("InChI=1S/C5000", "resource_limit", "numeric_limit"),
        ("InChI=1S/9999H2O", "resource_limit", "numeric_limit"),
        ("InChI=1S/2000H2O", "resource_limit", "formula_atom_limit"),
        ("InChI=1S/C/H1-" + "9" * 100, "resource_limit", "numeric_limit"),
        ("InChI=1S/" + "C" * 8192, "resource_limit", "input_byte_limit"),
    ],
)
def test_per_row_preflight_refuses_before_native(text, status, code, monkeypatch):
    def forbidden(*args):
        pytest.fail("rejected input must not reach native parser")

    monkeypatch.setattr(implementation, "from_inchi", forbidden)
    (row,) = ch.inchi_keys_checked([text])
    assert (row.index, row.source_text, row.status, row.code, row.key) == (
        0,
        text,
        status,
        code,
        None,
    )


def test_identity_difference_never_reads_key(monkeypatch):
    class Changed:
        standard_inchi = "InChI=1S/CH4/h1H4"

        @property
        def standard_inchikey(self):
            pytest.fail("must not return a key after changed readback")

    monkeypatch.setattr(implementation, "from_inchi", lambda text: Changed())
    (row,) = ch.inchi_keys_checked([WATER])
    assert (row.status, row.code, row.key) == ("unsupported", "identity_changed", None)


def test_native_failures_do_not_drop_rows(monkeypatch):
    parser = implementation.from_inchi

    def parse(text):
        if text == ETHANOL:
            raise ch.ChematicInputError("declared parse failure")
        return parser(text)

    monkeypatch.setattr(implementation, "from_inchi", parse)
    rows = ch.inchi_keys_checked([WATER, ETHANOL, WATER])
    assert [row.status for row in rows] == ["ok", "invalid", "ok"]
    assert rows[1].code == "inchi_parse" and rows[1].key is None


def test_backend_failure_and_invalid_key_are_explicit(monkeypatch):
    class NoBackend:
        @property
        def standard_inchi(self):
            raise RuntimeError("native backend unavailable")

    monkeypatch.setattr(implementation, "from_inchi", lambda text: NoBackend())
    (row,) = ch.inchi_keys_checked([WATER])
    assert (row.status, row.code) == ("conversion_error", "native_inchi_error")

    class BadKey:
        standard_inchi = WATER
        standard_inchikey = "not-a-key"

    monkeypatch.setattr(implementation, "from_inchi", lambda text: BadKey())
    (row,) = ch.inchi_keys_checked([WATER])
    assert (row.status, row.code, row.key) == ("conversion_error", "invalid_key", None)


@pytest.mark.parametrize(
    "error", [MemoryError("oom"), OSError("io"), KeyboardInterrupt()]
)
def test_unexpected_errors_propagate(error, monkeypatch):
    def fail(text):
        raise error

    monkeypatch.setattr(implementation, "from_inchi", fail)
    with pytest.raises(type(error)) as caught:
        ch.inchi_keys_checked([WATER])
    assert caught.value is error


@pytest.mark.parametrize(
    "values", [[WATER, None], [WATER, b"bytes"], WATER, b"bytes", None]
)
def test_bad_batch_types_precede_conversion(values, monkeypatch):
    def forbidden(*args):
        pytest.fail("batch validation must finish before native conversion")

    monkeypatch.setattr(implementation, "from_inchi", forbidden)
    with pytest.raises(TypeError):
        ch.inchi_keys_checked(values)


def test_iterator_error_precedes_conversion(monkeypatch):
    error = LookupError("iterator failed")

    def values():
        yield WATER
        raise error

    def forbidden(*args):
        pytest.fail("native work must not start")

    monkeypatch.setattr(implementation, "from_inchi", forbidden)
    with pytest.raises(LookupError) as caught:
        ch.inchi_keys_checked(values())
    assert caught.value is error


@pytest.mark.parametrize(
    "limit,error",
    [
        (True, TypeError),
        (np.bool_(True), TypeError),
        (1.0, TypeError),
        (0, ValueError),
        (10001, ValueError),
    ],
)
def test_limit_validation_before_iteration(limit, error):
    def forbidden():
        pytest.fail("limit validation must precede source iteration")
        yield WATER

    with pytest.raises(error):
        ch.inchi_keys_checked(forbidden(), max_records=limit)


def test_count_and_aggregate_caps_precede_conversion(monkeypatch):
    consumed = []

    def source():
        for i in range(10):
            consumed.append(i)
            yield WATER

    def forbidden(*args):
        pytest.fail("native work must not start")

    monkeypatch.setattr(implementation, "from_inchi", forbidden)
    with pytest.raises(ValueError):
        ch.inchi_keys_checked(source(), max_records=2)
    assert consumed == [0, 1, 2]
    with pytest.raises(ValueError, match="16 MiB"):
        ch.inchi_keys_checked(["x" * (8 << 20), "é" * (5 << 20)])


def test_empty_input_is_empty_and_numpy_integer_limit_is_valid():
    assert ch.inchi_keys_checked([], max_records=np.int64(1)) == ()


def test_native_error_message_is_bounded(monkeypatch):
    def fail(text):
        raise RuntimeError("x" * 5000)

    monkeypatch.setattr(implementation, "from_inchi", fail)
    (row,) = ch.inchi_keys_checked([WATER])
    assert row.message == "x" * 1024


def test_malformed_middle_identifier_does_not_hide_later_success():
    rows = ch.inchi_keys_checked([WATER, "InChI=1S/not-a-formula", ETHANOL])
    assert [row.status for row in rows] == ["ok", "invalid", "ok"]
    assert rows[2].key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


@pytest.mark.parametrize(
    "category,status",
    [
        ("resource_limit", "resource_limit"),
        ("unsupported", "unsupported"),
        ("ambiguous", "unsupported"),
    ],
)
def test_native_input_categories_are_not_mislabeled(category, status, monkeypatch):
    error = ch.ChematicInputError("declined")
    error.category = category
    error.code = "specific_native_reason"

    def fail(text):
        raise error

    monkeypatch.setattr(implementation, "from_inchi", fail)
    (row,) = ch.inchi_keys_checked([WATER])
    assert (row.status, row.code) == (status, "specific_native_reason")


@pytest.mark.parametrize("formula", ["C0", "0H2O", "H" + "0" * 5000])
def test_zero_and_long_zero_formula_fields_stay_record_errors(formula):
    (row,) = ch.inchi_keys_checked(["InChI=1S/" + formula])
    assert (row.status, row.code, row.key) == ("invalid", "invalid_formula", None)
