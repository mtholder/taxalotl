#!/usr/bin/env python3

from sqlalchemy import (Boolean,
                        Column,
                        ForeignKey,
                        Integer,
                        String,
                        Text,
                        )

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class Dataset(Base):
    __tablename__ = 'dataset'
    pk = Column(Integer, primary_key=True)
    name = Column(String)


class Authority(Base):
    __tablename__ = 'authority'
    pk = Column(Integer, primary_key=True)
    label = Column(String)


class NameFragment(Base):
    """A single word that can be a taxon name or part of a name"""
    __tablename__ = 'name_frag'
    pk = Column(Integer, primary_key=True)
    word = Column(String)


class TaxonName(Base):
    __tablename__ = 'taxon_name'
    pk = Column(Integer, primary_key=True)
    word_1_fk = Column(Integer, ForeignKey('name_frag.pk'), nullable=False)
    word_2_fk = Column(Integer, ForeignKey('name_frag.pk'), nullable=True)
    word_3_fk = Column(Integer, ForeignKey('name_frag.pk'), nullable=True)
    word_4_fk = Column(Integer, ForeignKey('name_frag.pk'), nullable=True)
    authority_fk = Column(Integer, ForeignKey('authority.pk'), nullable=True)


class TaxUUID(Base):
    __tablename__ = 'taxon_uuid'
    pk = Column(Integer, primary_key=True)
    uuid = Column(String, unique=True)


class Modified(Base):
    __tablename__ = 'modtime'
    pk = Column(Integer, primary_key=True)
    modified = Column(String, unique=True)


class AccordingTo(Base):
    __tablename__ = 'according_to'
    pk = Column(Integer, primary_key=True)
    according_to = Column(String, unique=True)


class Description(Base):
    __tablename__ = 'description'
    pk = Column(Integer, primary_key=True)
    description = Column(Text, unique=True)


class TaxonConcept(Base):
    __tablename__ = 'taxon_concept'
    pk = Column(Integer, primary_key=True)
    concept_id = Column(String, unique=True)


class SciNameUID(Base):
    __tablename__ = 'sci_name_uid'
    pk = Column(Integer, primary_key=True)
    sci_name_uid = Column(String, unique=True)


class CoLLink(Base):
    __tablename__ = 'col_link'
    pk = Column(Integer, primary_key=True)
    prim_uuid = Column(Integer, ForeignKey('taxon_uuid.pk'), nullable=True)
    syn_uuid = Column(Integer, ForeignKey('taxon_uuid.pk'), nullable=True)


class Rank(Base):
    __tablename__ = 'rank'
    pk = Column(Integer, primary_key=True)
    rank = Column(String)
    sorting_number = Column(Integer, nullable=True)


class TaxonSuppInfo(Base):
    __tablename__ = 'taxon_supp_info'
    pk = Column(Integer, primary_key=True)
    modified_fk = Column(Integer, ForeignKey('modtime.pk'), nullable=True)
    according_to_fk = Column(Integer, ForeignKey('according_to.pk'), nullable=True)
    description_fk = Column(Integer, ForeignKey('description.pk'), nullable=True)
    taxon_concept_fk = Column(Integer, ForeignKey('taxon_concept.pk'), nullable=True)
    col_link_fk = Column(Integer, ForeignKey('col_link.pk'), nullable=True)
    dataset_fk = Column(Integer, ForeignKey('dataset.pk'), nullable=True)
    uuid_fk = Column(Integer, ForeignKey('taxon_uuid.pk'), nullable=True)
    sci_name_uid_fk = Column(Integer, ForeignKey('sci_name_uid.pk'), nullable=True)
    name_status_fk = Column(Integer, ForeignKey('valid_name_status.pk'), nullable=True)
    syn_name_status_fk = Column(Integer, ForeignKey('syn_status.pk'), nullable=True)
    is_extinct = Column(Boolean)


class SynonymStatus(Base):
    __tablename__ = 'syn_status'
    pk = Column(Integer, primary_key=True)
    label = Column(String)


class ValidNameStatus(Base):
    __tablename__ = 'valid_name_status'
    pk = Column(Integer, primary_key=True)
    label = Column(String)


class Taxon(Base):
    __tablename__ = 'taxon'
    pk = Column(Integer, primary_key=True)
    name_fk = Column(Integer, ForeignKey('taxon_name.pk'))
    rank_fk = Column(Integer, ForeignKey('rank.pk'), nullable=True)
    tax_sup_info_fk = Column(Integer, ForeignKey('taxon_supp_info.pk'), nullable=True)


class Synonym(Base):
    __tablename__ = 'synonym'
    pk = Column(Integer, primary_key=True)
    name_fk = Column(Integer, ForeignKey('taxon_name.pk'))
    valid_taxon_fk = Column(Integer, ForeignKey('taxon.pk'))
    tax_sup_info_fk = Column(Integer, ForeignKey('taxon_supp_info.pk'), nullable=True)


class Edge(Base):
    __tablename__ = 'edge'
    pk = Column(Integer, primary_key=True)
    parent_fk = Column(Integer, ForeignKey('taxon.pk'))
    child_fk = Column(Integer, ForeignKey('taxon.pk'))
