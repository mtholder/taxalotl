#!/usr/bin/env python

from xml.sax import make_parser
from xml.sax.saxutils import XMLFilterBase
import logging

_LOG = logging.getLogger(__name__)


class WikispeciesPagesMetaCurrentParser(XMLFilterBase):
    """SAX parser for a pages-meta-current.xml dump from Wikispecies"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.in_page = False
        self.num_elements = 0
        self.current_tag = None
        self.current_title = None
        self.current_text = None

    def startElement(self, name, attrs):
        self.current_tag = name
        # print(f"{self.in_page} {name}")
        if name == "page":
            self.in_page = True
            self.current_title = None
            self.current_text = None
        self.num_elements += 1

    def endElement(self, name):
        if not self.in_page:
            return
        if name == "page":
            self.in_page = False
        elif name == "text":
            if self.current_text is None:
                return
            # if (
            #     ("{{Taxonbar" not in self.current_text)
            #     and ("Taxonavigation}}" not in self.current_text)
            #     and ("#REDIRECT [[" not in self.current_text.upper())
            # ):
            #     _LOG.debug(f"Non-taxon page {self.current_title.strip()}")

    def characters(self, content):
        if not self.in_page:
            return
        if self.current_tag == "title":
            if self.current_title is None:
                self.current_title = content
            else:
                self.current_title += content
        elif self.current_tag == "text":
            if self.current_text is None:
                self.current_text = content
            else:
                self.current_text += content


def parse_wikispecies(infp):
    reader = WikispeciesPagesMetaCurrentParser(make_parser())
    reader.parse(infp)
