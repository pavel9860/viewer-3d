"""Exception hierarchy for the CATIA V5 decoder."""

from __future__ import annotations


class CatV5Error(Exception):
    """Base class for every decoder error."""


class NotCatiaFile(CatV5Error):
    """The input does not start with the ``V5_CFV2`` container magic."""


class FormatError(CatV5Error):
    """A structure that was positively identified does not decode under its grammar.

    ``offset`` is the absolute file offset where decoding stopped.
    """

    def __init__(self, message: str, offset: int | None = None) -> None:
        self.offset = offset
        super().__init__(message if offset is None else f"{message} (at offset {offset})")
