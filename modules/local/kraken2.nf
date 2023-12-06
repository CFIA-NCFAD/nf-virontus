process KRAKEN2_PREPAREINDEX {
  input:
  path(db)

  output:
  path("kraken2-index/")

  script:
  """
  if [ -d "$db" ]; then
    if [[ (! -f "$db/hash.k2d") || (! -f "$db/opts.k2d") || (! -f "$db/taxo.k2d")]]; then
      echo "Missing Kraken 2 index files 'hash.k2d', 'opts.k2d' and/or 'taxo.k2d' in '$db'"
      exit 1
    fi
    ln -s $db kraken2-index
  elif [[ ("$db" =~ ^.*\\.tar\\.gz\$) && (\$(file -bL --mime-type "$db") = "application/gzip") ]]; then
    mkdir -p kraken2-index
    tar -xvzf $db -C kraken2-index
  else
    echo "Kraken 2 index '$db' is not a directory or a .tar.gz file."
    exit 1
  fi
  """
}
