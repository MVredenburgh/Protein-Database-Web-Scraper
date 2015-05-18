#!/usr/bin/env perl -w

use SequenceIterator qw(iterseq);
use Term::ANSIColor;
use Cwd;
use Carp;

my $rnafold = "RNAfold -T 21";

my %aa = ( 'ttt'=>'F',  'tct'=>'S',  'tat'=>'Y',  'tgt'=>'C',
	   'ttc'=>'F',  'tcc'=>'S',  'tac'=>'Y',  'tgc'=>'C',
	   'tta'=>'L',  'tca'=>'S',  'taa'=>'!',  'tga'=>'!',
	   'ttg'=>'L',  'tcg'=>'S',  'tag'=>'!',  'tgg'=>'W',
	   
	   'ctt'=>'L',  'cct'=>'P',  'cat'=>'H',  'cgt'=>'R',
	   'ctc'=>'L',  'ccc'=>'P',  'cac'=>'H',  'cgc'=>'R',
	   'cta'=>'L',  'cca'=>'P',  'caa'=>'Q',  'cga'=>'R',
	   'ctg'=>'L',  'ccg'=>'P',  'cag'=>'Q',  'cgg'=>'R',
	   
	   'att'=>'I',  'act'=>'T',  'aat'=>'N',  'agt'=>'S',
	   'atc'=>'I',  'acc'=>'T',  'aac'=>'N',  'agc'=>'S',
	   'ata'=>'I',  'aca'=>'T',  'aaa'=>'K',  'aga'=>'R',
	   'atg'=>'M',  'acg'=>'T',  'aag'=>'K',  'agg'=>'R',
	   
	   'gtt'=>'V',  'gct'=>'A',  'gat'=>'D',  'ggt'=>'G',
	   'gtc'=>'V',  'gcc'=>'A',  'gac'=>'D',  'ggc'=>'G',
	   'gta'=>'V',  'gca'=>'A',  'gaa'=>'E',  'gga'=>'G',
	   'gtg'=>'V',  'gcg'=>'A',  'gag'=>'E',  'ggg'=>'G' );


die "Usage: $0 input_dir output_dir" unless @ARGV == 2;
evaldir (@ARGV);

sub readFASTA {
    my ($filename) = @_;
    my @fasta;
    iterseq ($filename,
	     sub {
		 my ($name, $seq) = @_;
		 push @fasta, $seq;
	     });
    return @fasta;
}

sub readSingleFASTA {
    my ($filename) = @_;
    my @fasta = readFASTA ($filename);
    return $fasta[0];
}

sub fold {
    my ($seq) = @_;
    my $out = `echo $seq | $rnafold | sed -n '2p'`;
    chomp $out;
    $out =~ s/ .*//;
    $out =~ tr/\(\)/<>/;
    my $fold = "";
    for (my $i = 0; $i < length($out); ++$i) {
	my $c = substr($out,$i,1);
	$fold .= $c eq '.' ? substr($seq,$i,1) : $c;
    }
    return $fold;
}

sub getTATApos
{
    my ($dna) = @_;
    if ($dna =~ /tataat([acgt]{3,}aggagg[acgt]{6,10})(ATG|GTG|TTG)/g) {
	return pos($dna) - length($1) - length($2) - 7;
    } elsif ($dna =~ /tataat([acgt]{15,})(ATG|GTG|TTG)/g) {
	return pos($dna) - length($1) - length($2) - 7;
    } elsif ($dna =~ /tataat/g) {
	return pos($dna) - 7;
    } else {
	return undef;
    }
}

sub get35pos
{
    my ($dna) = @_;
    if ($dna =~ /ttgaca([acgt]{19}tataat)/g) {
	return pos($dna) - length($1) - 6;
    } elsif ($dna =~ /ttgaca/g) {
	return pos($dna) - 6;
    } else {
	return undef;
    }
}



sub getShineDalgarnoIndex
{
    my ($dna) = @_;
    if ($dna =~ /aggagg([acgt]{6,10})(ATG|GTG|TTG)/g) {
	return pos($dna) - 6 - length($1) - length($2);
    } elsif ($dna =~ /aggagg/g) {
	return pos($dna) - 6;
    } else {
	return undef;
    }
}

sub proportionShineDalgarnoExposed {
    my ($dna, $index) = @_;
    my $ss = fold($dna);
    print color('yellow'),substr(lc($ss),0,$index-1),colored['black on_yellow'],substr(uc($ss),$index,6),colored['yellow on_black'],substr(lc($ss),$index+6),color('reset'),"\n";
    my $sd = substr (lc($ss), $index, 6);
    my $n = 0;
    for (my $i = 0; $i < 6; ++$i) { ++$n if substr($sd,$i,1) =~ /[acgu]/ }
    return $n / 6;
}

sub maskEffector {
    my ($dna, $effector) = @_;
    my $x = "x" x length($effector);
    $dna = lc($dna);
    my $effector_rev = revcomp(lc($effector));
    $dna =~ s/$effector_rev/$x/g;
    return $dna;
}


sub getStartIndex
{
    my ($dna) = @_;
    if ($dna =~ /(ATG|GTG|TTG)/g) {
	return pos($dna) - 3;
    } else {
	return undef;
    }
}


sub getStopIndex
{
    my ($dna) = @_;
    if ($dna =~ /(TAG|TAA|TGA)/g) {
	return pos($dna) - 3;
    } else {
	return 0;
    }
}



sub getTermPos
{
    my ($dna) = @_;

    if ($dna =~ /(TAG|TAA|TGA)[acgt]+(t{7,})/g) {
	return pos($dna) - length($2);
    } elsif ($dna =~ /(t{7,})/) {
	return pos($dna) - length($1);
    } else {
	return undef;
    }

}


# Check restriction sites
sub testNoRestrictionSites {
    my ($seq, $sitesArrayRef) = @_;
    for my $site (@$sitesArrayRef) {
	$site = lc $site;
	$site = tr/u/t/;
	my $rev = revcomp ($site);
	if ($seq =~ /$site/ || $seq =~ /$rev/) {
	    return 0;
	}
    }
    return 1;
}


# Reverse complement
sub revcomp {
    my ($seq) = @_;
    my $rev = reverse $seq;
    $rev =~ tr/acgtuACGTU/tgcaaTGCAA/;
    return $rev;
}


# mean function
sub mean {
    my @x = @_;
    my ($x_sum, $x_n) = (0, 0);
    for my $x (@x) { $x_sum += $x; ++$x_n }
    return $x_n ? ($x_sum / $x_n) : 0;
}

# min function
sub min {
    my @x = @_;
    my $x_min = shift @x;
    for my $x (@x) { $x_min = $x if $x < $x_min }
    return $x_min;
}

# max function
sub max {
    my @x = @_;
    my $x_max = shift @x;
    for my $x (@x) { $x_max = $x if $x > $x_max }
    return $x_max;
}


# log score

sub score {
    my ($score, $test, $pass, $fail) = @_;
    carp "Undefined \$pass" unless defined($pass);
    $fail = $pass unless defined $fail;
    my ($points, $ret);
    $ret = $test == 0 ? 0 : ($test > 0 && $test < 1 ? int($test*$score + .5) : $score);
    if ($score < 0) {
	$ret = $score - $ret;
	$points = -$score . " points deducted";
    } else {  # $score > 0
	$points = $score . " points awarded";
    }
    if ($test < 0 || $test >= 1) {
	print colored ['bold white on_green'], "PASSED", color('bold yellow'), " ($ret/$points)", colored ['bold white on_green'], ": ", $pass, color('reset'), "\n";
    } elsif ($test > 0 || $test < 1) {
	print colored ['bold black on_cyan'], "PARTIAL", color('bold yellow'), " ($ret/$points)", colored ['bold black on_cyan'], ": ", $pass, color('reset'), "\n";
    } else {
	print colored ['white on_red'], "FAILED", color('bold yellow'), " ($ret/$points)", colored ['white on_red'], ": ", $fail, color('reset'), "\n";
    }
    return $ret;
}

# main function

sub evaldir {
    my ($srcdir, $dir) = @_;

    # symlink input files to output dir
    for my $filename (qw(protein.fasta sites.fasta effector.fasta codonfreq.txt params.txt)) {
	link "$srcdir/$filename", "$dir/$filename";
    }

    # read inputs & do basic preprocessing
    # protein.fasta
    my @protein = map (lc, readFASTA ("$srcdir/protein.fasta"));
    my $protlen = length $protein[0];
    grep (length != $protlen, @protein) && warn "Not all the proteins are the same length!";
    my $n_proteins = @protein;

    # sites.fasta
    my @sites = map (lc, readFASTA ("$srcdir/sites.fasta"));

    # codonfreq.txt
    my (%target_freq, %target_freq_norm);
    local *CODON;
    open CODON, "<$srcdir/codonfreq.txt";
    while (<CODON>) {
	my ($codon, $freq) = split;
	$target_freq{$codon} = $freq;
	$target_freq_norm{$aa{lc $codon}} += $freq;
    }
    close CODON;
    for my $codon (keys %target_freq) { $target_freq{$codon} /= $target_freq_norm{$aa{lc $codon}} }

    # params.txt
    local *PARAM;
    open PARAM, "<$srcdir/params.txt";
    my $N = <PARAM>;
    $N += 0;
    close PARAM;

    # unzip tarball & run project
    -e "$dir/project.tar.gz" && system "cd $dir; tar xvzf project.tar.gz";
    -e "$dir/runproject" && system "cd $dir; ./runproject";

    # read outputs
    # owneffector.fasta?
    my $own_effector = -e "$dir/owneffector.fasta";
    my $effector = readSingleFASTA ($own_effector ? "$dir/owneffector.fasta" : "$srcdir/effector.fasta");
    # output.fasta
    my @output = readFASTA ("$dir/output.fasta");

    # initialize scores
    my ($score, $rankScore) = (0, 0);
    $score += score (-5, !$own_effector, "Using contest-supplied effector.fasta", "Using team's owneffector.fasta");

    # Loop over the output sequences
    my $total = 0;
    my (@trans, @deltaFreq, @codonSep);
    print colored ['white on_blue'], "Looping over $N sequences", color('reset'), "\n";
    for (my $i = 0; $i < $N; ++$i) {
	my $seq = $output[$i];
	print color('magenta'), "Checking sequence #", $i+1, ".", color('reset'), "\n";
	print colored ['cyan on_blue'], $seq, color('reset'), "\n";
	my $is_rna = $seq =~ /[uU]/ ? 1 : 0;
	$seq =~ tr/uU/tT/;

	# Start & stop codon positions are flagged in UPPER CASE, so get them first, to constrain search.
	my $startpos = getStartIndex($seq);
	my $stoppos = getStopIndex($seq);

	# Translate, keeping track of codon usage & identical codon separation
	my $truncated = 0;
	my $prot = "";
	my (%lastPos, %codonFreq, %codonFreqNorm);
	my ($sepSum, $sepCount) = (0, 0);
	if (defined($startpos) && defined($stoppos)) {
	    for (my $i = $startpos; $i < $stoppos; $i += 3) {
		my $cpos = ($i - $startpos) / 3;
		my $codon = lc (substr ($seq, $i, 3));
		my $aa = lc $aa{$codon};
		if ($aa eq "!") { $truncated = 1 }
		else { $prot .= $aa }
		++$codonFreq{$codon};
		++$codonFreqNorm{$aa{$codon}};
		if (exists $lastPos{$codon}) {
		    $sepSum += $cpos - $lastPos{$codon};
		    ++$sepCount;
		}
		$lastPos{$codon} = $cpos;
	    }
	    my ($delta_freq_sum, $delta_freq_n) = (0, 0);
	    for my $codon (keys %codonFreq) {
		my $f = $codonFreq{$codon} / $codonFreqNorm{$aa{$codon}};
		$delta_freq_sum += $f - $target_freq{$codon};
	    }
	    push @trans, $prot;
	    push @deltaFreq, $delta_freq_sum / keys(%codonFreq);
	    push @codonSep, $sepCount ? ($sepSum / $sepCount) : 0;
	}

	# Check for Start Codon (5pts).
	print color('magenta'), "Checking for start codon (5 points).", color('reset'), "\n";
	$total += score (1, length($seq) > 0, "Found sequence #".($i+1));
	$total += score (1, $seq =~ /^[acgt]+$/i ? 1 : 0, "Looks like nucleotide sequence");
	$total += score (1, !$is_rna, "No RNA U's");
	$total += score (1, $seq =~ /^[acgt]*[ACGT]{3}[acgt]*[ACGT]{3}[acgt]*$/ ? 1 : 0, "Exactly two codons are annotated");
	$total += score (1, defined($startpos), "Found start codon");

	# Check for Stop Codon (5pts).
	print color('magenta'), "Checking for stop codon (5 points).", color('reset'), "\n";
	$total += score (1, defined($stoppos), "Found stop codon");
	$total += score (1, defined($startpos) && defined($stoppos) && $stoppos > $startpos, "Stop codon is after start codon");
	$total += score (1, defined($startpos) && defined($stoppos) && $stoppos > $startpos && (($stoppos - $startpos) % 3) == 0, "Separation between start & stop codons is multiple of 3 nucleotides");
	$total += score (1, defined($startpos) && defined($stoppos) && $stoppos > $startpos && $stoppos == $startpos + $protlen * 3, "Translated protein has correct length");
	$total += score (1, defined($startpos) && defined($stoppos) && $stoppos > $startpos && (($stoppos - $startpos) % 3) == 0 && !$truncated, "No premature stop codons");

	# Check for Shine-Dalgarno sequence (5pts).
	print color('magenta'), "Checking for Shine-Dalgarno sequence (5 points).", color('reset'), "\n";
	my $sdpos = getShineDalgarnoIndex($seq);
	$total += score (2, defined($sdpos), "Found consensus Shine-Dalgarno sequence");
	$total += score (1, defined($sdpos) && defined($startpos) && $startpos > $sdpos, "Shine-Dalgarno sequence is before start codon");
	$total += score (1, defined($sdpos) && defined($startpos) && $sdpos+6 >= $startpos - 10 && $sdpos+6 <= $startpos - 6, "Shine-Dalgarno sequence is 6-10 bases upstream of start codon", "Shine-Dalgarno sequence is ".($startpos-$sdpos-6)." bases upstream of start codon (expected 6-10)");
	$total += score (1, defined($sdpos) && defined($startpos) && !(substr($seq,$sdpos,$startpos-$sdpos) =~ /(atg|gtg|ttg)/), "Annotated start codon is first to follow Shine-Dalgarno sequence, consistent with 'ribosome scanning'");

	# Check for Promoter region (5pts).
	print color('magenta'), "Checking for basic bacterial promoter (5 points).", color('reset'), "\n";
	my $tatapos = getTATApos($seq);
	my $pos35 = get35pos($seq);
	$total += score (2, defined($tatapos), "Found consensus Pribnow box");
	$total += score (1, defined($tatapos) && defined($sdpos) && $sdpos - $tatapos >= 10, "Pribnow box is at least 10 bases upstream of Shine-Dalgarno site");
	$total += score (1, defined($pos35), "Found consensus -35 box");
	$total += score (1, defined($tatapos) && defined($pos35) && ($tatapos - $pos35 >= 23 && $tatapos - $pos35 < 27), "Pribnow box is 23-27 bases downstream of -35 box");

	# Check for Terminator (5pts).
	print color('magenta'), "Checking for transcriptional terminator (5 points).", color('reset'), "\n";
	my $termpos = getTermPos($seq);
	my $termss = defined($termpos) ? fold(substr($seq,0,$termpos)) : "";
	$total += score (1, defined($termpos), "Found poly-T tail sequence of transcription terminator");
	$total += score (1, defined($termpos) && $termpos >= $stoppos + 3, "Poly-T terminator sequence is downstream of stop codon");
	if (defined $termpos) {
	    print color('yellow'),substr(lc($termss),0,$termpos-1),colored['black on_yellow'],substr(uc($seq),$termpos,7),colored['yellow on_black'],substr(lc($seq),$termpos+7),color('reset'),"\n";
	}
	$total += score (3, $termss =~ /<{4,}[acgt]+>{4,}[acgt]*$/ ? 1 : 0, "Found terminator stem-loop");

	# Check for functional Riboswitch (20pts: half if Shine-Dalgarno sequence exposed when effector present, half if exposed when effector absent).
	print color('magenta'), "Checking for functional riboswitch (20 points).", color('reset'), "\n";
	if ($sdpos > 0) {
	    $total += score (10, proportionShineDalgarnoExposed($seq,$sdpos) < 1, "Shine-Dalgarno occluded in secondary structure when effector absent");
	    my $maskedSeq = maskEffector($seq,$effector);
	    my $efflen = length($effector);
	    my $gotMask = $maskedSeq =~ /x{$efflen}/;
	    $total += score (4, $gotMask, "Found exact reverse-complement of effector in sequence");
	    $total += score (6, $gotMask ? proportionShineDalgarnoExposed($maskedSeq,$sdpos) : 0, "Shine-Dalgarno exposed in secondary structure when effector present");
	} else {
	    score (20, 0, undef, "Couldn't find Shine-Dalgarno sequence, so no riboswitch points");
	}

	# Check for Restriction sites and their reverse complement (10pts).
	print color('magenta'), "Checking for restriction sites (10 points).", color('reset'), "\n";
	$total += score (10, testNoRestrictionSites($seq), "No restriction sites");
    }

    print colored ['white on_blue'], "End loop over $N sequences", color('reset'), "\n";
    print "\n", colored ['white on_blue'], "Base score + (feature score averaged over $N sequences): ", color('cyan'), $score, " + (", $total, " / ", $N, ") = ", $score + $total / $N, color('reset'), "\n";
    $score += $total / $N;

    # Perform competitive ranking against other teams (up to 20pts depending on rank).
    # RankScore is a sum of the following:
    # S / (sequence length), the length-normalized total column entropy in an alignment of all your encoded protein sequences (rewarding diversity in your generated protein sequences).
    my $S_sum = 0;
    for (my $col = 0; $col < $protlen; ++$col) {
	my %dist;
	my $n = 0;
	for my $trans (@trans) {
	    if (length($trans) > $col) {
		my $aa = substr($trans,$col,1);
		++$dist{$aa};
		++$n;
	    }
	}
	for my $aa_count (values %dist) {
	    my $p = $aa_count / $n;
	    if ($p > 0) {
		$S_sum -= $p * log($p);
	    }
	}
    }
    my $S = $S_sum / $protlen;

    # -M, a penalty term based on the maximum percentage identity between any one of your proteins and any protein from protein.fasta (rewarding generation of new sequences, rather than exact copies of the sequence in that file).
    my $M = 0;
    for my $trans (@trans) {
	for my $prot (@protein) {
	    my ($id_sum, $id_count) = (0, 0);
	    for (my $i = 0; $i < $protlen; ++$i) {
		++$id_sum if lc(substr($trans,$i,1)) eq lc(substr($prot,$i,1));
		++$id_count;
	    }
	    my $id = $id_sum / $id_count;
	    $M = $id if $id > $M;
	}
    }

    # L / (sequence length) / (number of columns), the length- and depth-normalized total log-likelihood of all your protein sequences using a position-specific weight matrix model based on the alignment in protein.fasta (rewarding identification of the conserved/variable sites).
    my $L = 0;
    for (my $col = 0; $col < $protlen; ++$col) {
	my %dist;
	for my $prot (@protein) {
	    my $aa = substr($prot,$col,1);
	    ++$dist{$aa};
	}
	my %p = map (($_ => $dist{$_} / $n_proteins), keys %dist);
	for my $trans (@trans) {
	    if (length($trans) > $col) {
		my $aa = substr($trans,$col,1);
		$L += log $p{$aa};
	    } else {
		$L += log(.001);  # generic penalty
	    }
	}
    }
    $L /= $protlen;

    # -Delta_CodonFreq, the difference between the a codon's frequency in your output sequence and the target frequency, averaged across all codons (rewarding correct codon usage).
    my $deltaCodonFreq = mean (@deltaFreq);

    # -(Avg Space between similar codons) / (sequence length), the average space between identical/synonymous codons, averaged for all codons and divided by sequence length (rewarding avoidance of repeat codons).
    my $meanCodonSep = mean (@codonSep) / $protlen;

    # Print components
    print colored ['white on_blue'], sprintf("Rank score components: entropy=%.4f(+) maxID=%.4f(-) negLogLike=%.4f(-) deltaFreq=%.4f(-) codonSep=%.4f(-)",$S,$M,-$L,$deltaCodonFreq,$meanCodonSep), color('reset'), "\n";
    print colored ['white on_blue'], "...where (+) scores should be maximized, and (-) scores should be minimized", color('reset'), "\n";

    # Compute rank score
    $rankScore = sprintf ("%.4f", $S - $M + $L - $deltaCodonFreq - $meanCodonSep);

    local *SCORE;
    open SCORE, ">$dir/FeatureScore.txt";
    print SCORE $score, "\n";
    close SCORE;

    local *RANKSCORE;
    open RANKSCORE, ">$dir/RankScore.txt";
    print RANKSCORE $rankScore, "\n";
    close RANKSCORE;

    print "\n";
    print colored ['white on_blue'], "Feature score: ", color('cyan'), $score, color('reset'), "\n";
    print colored ['white on_blue'], "Rank score: ", color('cyan'), $rankScore, color('reset'), "\n";
}

