import React from "react";
import {
  createViewState,
  createJBrowseTheme,
  JBrowseLinearGenomeView,
  ThemeProvider,
} from "@jbrowse/react-linear-genome-view";

const theme = createJBrowseTheme();

const assembly = {
  name: "GRCh38",
  sequence: {
    type: "ReferenceSequenceTrack",
    trackId: "GRCh38-ReferenceSequenceTrack",
    adapter: {
      type: "BgzipFastaAdapter",
      fastaLocation: {
        uri:
          "http://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
      },
      faiLocation: {
        uri:
          "http://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai",
      },
      gziLocation: {
        uri:
          "http://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi",
      },
    },
  },
  aliases: ["hg38"],
  refNameAliases: {
    adapter: {
      type: "RefNameAliasAdapter",
      location: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/hg38_aliases.txt",
      },
    },
  },
};

const tracks = [
  {
    type: "FeatureTrack",
    trackId: "ncbi_refseq_109_hg38",
    name: "NCBI RefSeq (GFF3Tabix)",
    assemblyNames: ["GRCh38"],
    category: ["Annotation"],
    adapter: {
      type: "Gff3TabixAdapter",
      gffGzLocation: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/ncbi_refseq/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sorted.gff.gz",
      },
      index: {
        location: {
          uri:
            "https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/ncbi_refseq/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sorted.gff.gz.tbi",
        },
      },
    },
  },
  {
    type: "AlignmentsTrack",
    trackId: "shortread_cram",
    name: "shortread_cram",
    assemblyNames: ["GRCh38"],
    adapter: {
      type: "CramAdapter",
      cramLocation: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/shortread_test_snippet.cram",
      },
      craiLocation: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/shortread_test_snippet.cram.crai",
      },
      sequenceAdapter: {
        type: "BgzipFastaAdapter",
        fastaLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz",
        },
        faiLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.fai",
        },
        gziLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.gzi",
        },
      },
    },
  },
  {
    type: "AlignmentsTrack",
    trackId: "shortread_bam",
    name: "shortread_bam",
    assemblyNames: ["GRCh38"],
    adapter: {
      type: "BamAdapter",
      bamLocation: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/shortread_test_snippet.bam",
      },
      index: {
        location: {
          uri:
            "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/shortread_test_snippet.bam.bai",
        },
      },
      sequenceAdapter: {
        type: "BgzipFastaAdapter",
        fastaLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz",
        },
        faiLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.fai",
        },
        gziLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.gzi",
        },
      },
    },
  },
  {
    type: "AlignmentsTrack",
    trackId: "longread_cram",
    name: "longread_cram",
    assemblyNames: ["GRCh38"],
    adapter: {
      type: "CramAdapter",
      cramLocation: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/longread_test_snippet.cram",
      },
      craiLocation: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/longread_test_snippet.cram.crai",
      },
      sequenceAdapter: {
        type: "BgzipFastaAdapter",
        fastaLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz",
        },
        faiLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.fai",
        },
        gziLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.gzi",
        },
      },
    },
  },
  {
    type: "AlignmentsTrack",
    trackId: "longread_bam",
    name: "longread_bam",
    assemblyNames: ["GRCh38"],
    adapter: {
      type: "BamAdapter",
      bamLocation: {
        uri:
          "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/longread_test_snippet.bam",
      },
      index: {
        location: {
          uri:
            "https://s3.amazonaws.com/jbrowse.org/genomes/hg19/longread_test_snippet.bam.bai",
        },
      },
      sequenceAdapter: {
        type: "BgzipFastaAdapter",
        fastaLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz",
        },
        faiLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.fai",
        },
        gziLocation: {
          uri: "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz.gzi",
        },
      },
    },
  },
];

const defaultSession = {
  name: "My session",
  view: {
    id: "linearGenomeView",
    type: "LinearGenomeView",
    tracks: [
      {
        type: "ReferenceSequenceTrack",
        configuration: "GRCh38-ReferenceSequenceTrack",
        displays: [
          {
            type: "LinearReferenceSequenceDisplay",
            configuration:
              "GRCh38-ReferenceSequenceTrack-LinearReferenceSequenceDisplay",
          },
        ],
      },
      {
        type: "FeatureTrack",
        configuration: "ncbi_refseq_109_hg38",
        displays: [
          {
            type: "LinearBasicDisplay",
            configuration: "ncbi_refseq_109_hg38-LinearBasicDisplay",
          },
        ],
      },
    ],
  },
};

function View() {
  const state = createViewState({
    assembly,
    tracks,
    location: "10:29,838,737..29,838,819",
    defaultSession,
  });
  return (
    <ThemeProvider theme={theme}>
      <JBrowseLinearGenomeView viewState={state} />
    </ThemeProvider>
  );
}

export default View;
