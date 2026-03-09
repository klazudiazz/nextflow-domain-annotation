nextflow.enable.dsl = 2


// 1) extractorfs: z genomu wycina ORFy >=900bp i tłumaczy do białek (*.faa)

process extractorfs {
    conda "envs/powb-pyhmmer.yml"

    input:
        tuple val(asmid), path(genome)
        path alph
        path tab

    output:
        path "${asmid}.faa"

    script:
        """
        gunzip -c $genome | \
          extractorfs --alph ${alph} \
                      --tab ${tab} \
                      --asmacc ${asmid} \
                      --minlen 900 \
                      --trans "${asmid}.faa" \
          > /dev/null
        """
}

// 2) stage_faa: tworzy JEDEN katalog z linkami do wszystkich *.faa (1 task)

process stage_faa {
    conda "envs/powb-pyhmmer.yml"

    input:
        path faa_files

    output:
        path "faa_stage"

    script:
        """
        mkdir -p faa_stage

        for f in ${faa_files}; do
          cp -L "\$f" faa_stage/
        done
        """
}


// 3) uniquetrans: buduje niereundantny zbiór białek (unique.faa) + tsv z klastrami i długościami
process uniquetrans {
    publishDir "output/unique", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"

    input:
        path faa_dir

    output:
        path "unique.faa"
        path "unique_clusts.tsv"
        path "unique_lengths.tsv"

    script:
        """
        unique.py --input "${faa_dir}/*.faa" \
                  --output unique.faa
        """
}

// 4) splitfasta: dzieli unique.faa na N chunków
process splitfasta {
    publishDir "output/splitfasta", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"

    input:
        path uniq_faa
        val num_chunks

    output:
        path "chunk_*.faa"

    script:
        """
        splitfasta.py --chunks ${num_chunks} \
                      --input ${uniq_faa} \
                      --outdir .
        """
}

// 5) hmmfetch: wycina profil domeny PF01551 z Pfam-A.hmm.gz
process hmmfetch {
    publishDir "output/hmmfetch", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"

    input:
        path pfam_db
        val domain_acc

    output:
        path "domain.hmm"

    script:
        """
        hmmfetch.py --domacc ${domain_acc} \
                    --pfamdb ${pfam_db} \
                    --output domain.hmm
        """
}

// 6) hmmselsearch: szuka domeny M23 w chunkach unique.faa (po 1 TSV na chunk)
process hmmselsearch {
    publishDir "output/hmmselsearch", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"
    tag "${chunk.baseName}"

    input:
        path chunk
        path domain_hmm
        val domcov
        val iEval

    output:
        path "hmmselsearch_${chunk.baseName}.tsv"

    script:
        """
        hmmsearch.py --hmm ${domain_hmm} \
                     --seqs ${chunk} \
                     --domcov ${domcov} \
                     --iEval ${iEval} \
                     --output hmmselsearch_${chunk.baseName}.tsv
        """
}

// 7) extracttrans: wycina sekwencje białkowe zawierające M23 z ORFów (*.faa) po asmacc
process extracttrans {
    publishDir "output/extracttrans", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"
    tag "${hmm_results.baseName}"

    input:
        path hmm_results
        path seqdir

    output:
        path "extracttrans_${hmm_results.baseName}.faa"

    script:
        """
        extractfasta.py --seqdir ${seqdir} \
                        --hmmres ${hmm_results} \
                        --output extracttrans_${hmm_results.baseName}.faa
        """
}

// 8) hmmallsearch: szuka wszystkich domen PFAM w sekwencjach M23-pozytywnych
process hmmallsearch {
    publishDir "output/hmmallsearch", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"
    tag "${trans_faa.baseName}"

    input:
        path trans_faa
        path pfam_db
        val domcov
        val iEval

    output:
        path "hmmallsearch_${trans_faa.baseName}.tsv"

    script:
        """
        hmmsearch.py --hmm ${pfam_db} \
                     --seqs ${trans_faa} \
                     --domcov ${domcov} \
                     --iEval ${iEval} \
                     --output hmmallsearch_${trans_faa.baseName}.tsv
        """
}

// 9) annotdom: generuje GFF3 na podstawie zintegrowanych wyników PFAM i finalnego FASTA
process annotdom {
    publishDir "output/integrated", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"

    input:
        path allres_tsv
        path seqs_faa

    output:
        path "integrated.gff3"

    script:
        """
        annot.py --allres ${allres_tsv} \
                 --seqs ${seqs_faa} \
                 --output integrated.gff3
        """
}

// 10) archchart: renderuje HTML z architekturami domenowymi M23 (wymaga style.css i tmpl.html)
process archchart {
    publishDir "output/archchart", mode: 'rellink'
    conda "envs/powb-pyhmmer.yml"

    input:
        val domacc
        path style_file
        path tmpl_file
        path clstlen_tsv
        path hmmres_tsv

    output:
        path "archchart_${domacc}.html"

    script:
        """
        archchart.py --domacc ${domacc} \
                     --style ${style_file} \
                     --tmpl  ${tmpl_file} \
                     --clstlen ${clstlen_tsv} \
                     --hmmres ${hmmres_tsv} \
                     --output archchart_${domacc}.html
        """
}



workflow {

    // params.genom_seqs może być globem (np. "/path/*.fna.gz") albo katalogiem ("/path/")
    def genome_pattern = params.genom_seqs.toString().endsWith('/')
        ? "${params.genom_seqs}/*_genomic.fna*"
        : params.genom_seqs

    // kanał wejściowy genomów: tuple(asmid, plik)
    g_ch = Channel.fromPath(genome_pattern)
        .map { p ->
            def m = (p.getName() =~ /^(.+)_genomic\.fna(\.gz)?$/)
            if( !m ) error "Zła nazwa pliku: ${p.getName()}"
            tuple(m[0][1], p)
        }

    // 1) ORFy -> 118k plików *.faa (w work/)
    faa_ch = extractorfs(g_ch, file(params.alph), file(params.tab))

    // 2) stage: jeden katalog z linkami do wszystkich *.faa (żeby unique.py miał glob)
    faa_stage_dir = stage_faa( faa_ch.collect() )

    // 3) nieredundantny zbiór białek
    (uniq_faa_ch, uniq_clust_ch, uniq_len_ch) = uniquetrans(faa_stage_dir)

    // 4) split na chunki (pomija puste)
    chunks_ch = splitfasta(uniq_faa_ch, params.cpus)
        .flatten()
        .filter { it.size() > 0 }

    chunks_ch.ifEmpty { error "Brak niepustych chunków po splitfasta (unique.faa jest puste?)" }

    // 5) domena PF01551 wycięta z Pfam-A
    domain_hmm_ch = hmmfetch(file(params.pfam_db), params.domain_acc)

    // 6) wyszukuje M23 w chunkach
    hmmsel_ch = hmmselsearch(chunks_ch, domain_hmm_ch, params.domcov, params.iEval)

    // 7) wycina sekwencje M23-pozytywne z plików *.faa
    seqdir = faa_stage_dir

    extracttrans_ch = extracttrans(hmmsel_ch, seqdir)
        .filter { it.exists() && it.size() > 0 }

    // 8) wyszukuje wszystkie domeny PFAM w M23-pozytywnych sekwencjach
    hmmallsearch_ch = hmmallsearch(extracttrans_ch, file(params.pfam_db), params.domcov, params.iEval)

    // 9) integrate (2x collectFile) - pliki zbiorcze
    m23_all_faa = extracttrans_ch.collectFile(
        name: 'M23_sequences.faa',
        storeDir: 'output/integrate',
        keepHeader: true
    )

    pfam_all_tsv = hmmallsearch_ch.collectFile(
        name: 'PFAM_results.tsv',
        storeDir: 'output/integrate',
        keepHeader: true
    )

    // 10) GFF3 z anotacjami domen
    annotdom(pfam_all_tsv, m23_all_faa)

    // 11) HTML z architekturami domenowymi
    archchart(
        params.domain_acc,
        file(params.arch_style),
        file(params.arch_tmpl),
        uniq_len_ch,
        pfam_all_tsv
    )
}
