/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;

import com.google.errorprone.annotations.Var;

import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Combines multiple variant files into a single variant file.
 *
 * <h3>Inputs</h3>
 *  <ul>
 *      <li>One or more input file in VCF format (can be gzipped, i.e. ending in ".vcf.gz", or binary compressed, i.e. ending in ".bcf").</li>
 *      <li>Optionally a sequence dictionary file (typically name ending in .dict) if the input VCF does not contain a
 *          complete contig list and if the output index is to be created (true by default).</li>
 *  </ul>
 *  <p>
 *  The input variant data must adhere to the following rules:
 *     <ul>
 *         <li>If there are samples, those must be the same across all input files.</li>
 *         <li>Input file headers must be contain compatible declarations for common annotations (INFO, FORMAT fields) and filters.</li>
 *         <li>Input files variant records must be sorted by their contig and position following the sequence dictionary provided
 *         or the header contig list.</li>
 *     </ul>
 * </p>
 * <p>
 *     You can either directly specify the list of files by specifying <code>INPUT</code> multiple times, or provide a list
 *     in a file with name ending in ".list" to <code>INPUT</code>.
 * </p>
 *
 * <h3>Outputs</h3>
 * A VCF sorted (i) according to the dictionary and (ii) by coordiante. 
 * <h3>Usage examples</h3>
 * <h4>Example 1:</h4>
 *     We combine several variant files in different formats, where at least one of them contains the contig list in its header.
 * <pre>
 *     java -jar picard.jar MungeVcfs \
 *          I=input_variants.01.vcf \
 *          I=input_variants.02.vcf.gz \
 *          O=output_variants.vcf.gz
 * </pre>
 * <h4>Example 2:</h4>
 *      Similar to example 1 but we use an input list file to specify the input files:
 * <pre>
 *     java -jar picard.jar MungeVcfs \
 *          I=input_variant_files.list \
 *          O=output_variants.vcf.gz
 * </pre>
 *
 * @since 1.0.1
 */
@CommandLineProgramProperties(
		oneLineSummary = MungeVcfs.SUMMARY_FIRST_SENTENCE,
        summary = MungeVcfs.SUMMARY,
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature
public class MungeVcfs extends CommandLineProgram {

	static final String SUMMARY_FIRST_SENTENCE = "Combines multiple variant files into a single variant file";
	static final String SUMMARY = "<p>" + SUMMARY_FIRST_SENTENCE + ".</p>" + 
			"<h3>Inputs</h3>" + 
			"<ul>" + 
			"      <li>One or more input file in VCF format (can be gzipped, i.e. ending in \".vcf.gz\", or binary compressed, i.e. ending in \".bcf\").</li>" + 
			"      <li>Optionally a sequence dictionary file (typically name ending in .dict) if the input VCF does not contain a" + 
			"          complete contig list and if the output index is to be created (true by default).</li>" + 
			"  </ul>" + 
			"  <p>" + 
			"  The input variant data must adhere to the following rules:</p>" + 
			"     <ul>" + 
			"         <li>If there are samples, those must be the same across all input files.</li>" + 
			"         <li>Input file headers must be contain compatible declarations for common annotations (INFO, FORMAT fields) and filters.</li>" + 
			"         <li>Input files variant records must be sorted by their contig and position following the sequence dictionary provided" + 
			"         or the header contig list.</li>" + 
			"     </ul>" + 
			" <p>You can either directly specify the list of files by specifying <code>INPUT</code> multiple times, or provide a list" + 
			"     in a file with name ending in \".list\" to <code>INPUT</code>.</p>" + 
			" <h3>Outputs</h3>" + 
			" <p>A VCF sorted (i) according to the dictionary and (ii) by coordinate.</p>" +
			" <h3>Usage examples</h3>" + 
			" <h4>Example 1:</h4>" + 
			" <p>We combine several variant files in different formats, where at least one of them contains the contig list in its header.</p>" + 
			" <pre>java -jar picard.jar MungeVcfs \\\n" + 
			"          I=input_variants.01.vcf \\\n" + 
			"          I=input_variants.02.vcf.gz \\\n" + 
			"          O=output_variants.vcf.gz</pre>" + 
			" <h4>Example 2:</h4>" + 
			" <p>Similar to example 1 but we use an input list file to specify the input files:</p>" + 
			" <pre>java -jar picard.jar MungeVcfs \\\n" + 
			"          I=input_variant_files.list \\\n" + 
			"          O=output_variants.vcf.gz</pre><hr/>";  
	  	
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
              doc="VCF or BCF input files (File format is determined by file extension), or a file having a '.list' suffix containing the path to the files, one per line.", minElements=1)
    public List<PicardHtsPath> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The merged VCF or BCF file. File format is determined by file extension.")
    public File OUTPUT;

    @Argument(shortName = "D", doc = "The index sequence dictionary to use instead of the sequence dictionary in the input files", optional = true)
    public PicardHtsPath SEQUENCE_DICTIONARY;

    @Argument(doc = "Comment(s) to include in the merged output file's header.", optional = true, shortName = "CO")
    public List<String>  COMMENT = new ArrayList<>();

    private final static String SEQ_DICT_REQUIRED = "A sequence dictionary must be available (either through the input file or by setting it explicitly).";

    private final Log log = Log.getInstance(MergeVcfs.class);

    public MungeVcfs() {
        this.CREATE_INDEX = true;
    }

    @Override
    protected int doWork() {
        final ProgressLogger progress = new ProgressLogger(log, 10);
        // store sampleList in a set
        final Set<String> sampleSet = new HashSet<>();
        final List<Path> inputPaths = INPUT.stream().map(PicardHtsPath::toPath).collect(Collectors.toList());
        final List<Path> unrolledPaths = IOUtil.unrollPaths(inputPaths, FileExtensions.VCF_LIST.toArray(new String[]{}));
        final ArrayList<CloseableIterator<VariantContext>> iteratorCollection = new ArrayList<>(unrolledPaths.size());
        final Collection<VCFHeader> headers = new HashSet<>(unrolledPaths.size());
        VariantContextComparator variantContextComparator = null;
        SAMSequenceDictionary sequenceDictionary = null;

        if (SEQUENCE_DICTIONARY != null)
            sequenceDictionary = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(SEQUENCE_DICTIONARY.toPath()).getFileHeader().getSequenceDictionary();

        for (final Path path : unrolledPaths) {
            IOUtil.assertFileIsReadable(path);
            final VCFFileReader fileReader = new VCFFileReader(path, false);
            final VCFHeader fileHeader = fileReader.getFileHeader();
            if (fileHeader.getContigLines().isEmpty()) {
                if (sequenceDictionary == null) {
                    throw new IllegalArgumentException(SEQ_DICT_REQUIRED);
                } else {
                    fileHeader.setSequenceDictionary(sequenceDictionary);
                }
            }

            if (variantContextComparator == null) {
                variantContextComparator = fileHeader.getVCFRecordComparator();
            } else {
                if (!variantContextComparator.isCompatible(fileHeader.getContigLines())) {
                    throw new IllegalArgumentException(
                            "The contig entries in input path " + path.toAbsolutePath() + " are not compatible with the others.");
                }
            }

            if (sequenceDictionary == null) sequenceDictionary = fileHeader.getSequenceDictionary();
            
            List<String> samplesInCurrentFile = fileHeader.getSampleNamesInOrder();
            System.out.println("file: " + path.toAbsolutePath() + " has " + samplesInCurrentFile.size() + " samples");

            sampleSet.addAll(samplesInCurrentFile);
            
            // add comments in the first header
            if (headers.isEmpty()) {
                COMMENT.forEach(C -> fileHeader.addMetaDataLine(new VCFHeaderLine("MungeVcfs.comment", C)));
            }
            headers.add(fileHeader);
            iteratorCollection.add(fileReader.iterator());

        }
        System.out.println("sampleSet has " + sampleSet.size() + " samples");
        
        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new PicardException(String.format("Index creation failed. %s", SEQ_DICT_REQUIRED));
        }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary);

        if (CREATE_INDEX) {
            builder.setOption(Options.INDEX_ON_THE_FLY);
        } else {
            builder.unsetOption(Options.INDEX_ON_THE_FLY);
        }
        final VariantContextWriter writer = builder.build();

        writer.writeHeader(new VCFHeader(VCFUtils.smartMergeHeaders(headers, false), sampleSet));


        ArrayList<VariantContext> variantContexts = new ArrayList<>();
  
        for (final CloseableIterator<VariantContext> iterator : iteratorCollection) {
            if (iterator.hasNext()) {
                    VariantContext vc = iterator.next();
                    variantContexts.add(vc);
            }
        }
        //System.out.println("variantContexts has " + variantContexts.size() + " variant contexts");
        // while some variant contexts are not null
        Integer variantCount = 0;

        while (!allItemsAreNull(variantContexts)) {
            VariantContext firstVariantContext = findFirstVariantContext(variantContextComparator, variantContexts);

            //printVariantContexts(variantContexts);
            //System.out.println("firstVariantContext is " + firstVariantContext);


            VariantContextBuilder variantContextBuilder = new VariantContextBuilder(firstVariantContext);
            variantContextBuilder.noGenotypes();
            
            ArrayList<VariantContext> currentVariantContexts = new ArrayList<>();
            for (VariantContext vc : variantContexts) {
                if (variantContextComparator.compare(vc, firstVariantContext) == 0) {
                    currentVariantContexts.add(vc);
                }
            }
            //System.out.println("variantContexts has " + variantContexts.size() + " variant contexts");
            //System.out.println("currentVariantContexts has " + currentVariantContexts.size() + " variant contexts");

            List<Genotype> genotypes = new ArrayList<Genotype>(sampleSet.size());
            for (String sample : sampleSet)
            {
                // look through the variant contexts in the current position
                // and find the one that has the sample
                Genotype g = null;
                for (VariantContext vc : currentVariantContexts) {
                    g = vc.getGenotype(sample);

                    
                    if (g != null) {
                        break;
                    }
                }
                if (g == null)
                {
                    g = GenotypeBuilder.createMissing(sample, 2);
                }
                genotypes.add(g);
            }
            variantContextBuilder.genotypes(genotypes);
            VariantContext vc = variantContextBuilder.make();
            writer.add(vc);
            
            //variantContexts.remove(0);
            progress.record(firstVariantContext.getContig(), firstVariantContext.getStart());
            for (int i = 0; i < iteratorCollection.size(); i++) {
                CloseableIterator<VariantContext> iterator = iteratorCollection.get(i);
                VariantContext variantContext = variantContexts.get(i);
                if (variantContextComparator.compare(variantContext, firstVariantContext) <= 0) {
                    variantContexts.set(i, iterator.next());
                }
            }
            
            variantCount++;
            // if(variantCount > 10) {
            //     break;
            // }
        }

        writer.close();
        return 0;
    }

    private VariantContext findFirstVariantContext(VariantContextComparator variantContextComparator,
            ArrayList<VariantContext> variantContexts) {
        VariantContext firstVariantContext = null;
        for (VariantContext vc : variantContexts) {
            if (firstVariantContext == null) {
                firstVariantContext = vc;
            } else {
                if (variantContextComparator.compare(vc, firstVariantContext) < 0) {
                    firstVariantContext = vc;
                }
            }
        }
        return firstVariantContext;
    }

    private <T> Boolean allItemsAreNull(ArrayList<T> items) {
        for (T item : items) {
            if (item != null) {
                return false;
            }
        }
        return true;
    }

    private void printVariantContexts(ArrayList<VariantContext> variantContexts) {
        String variantContextsString = "";
        for (VariantContext vc : variantContexts) {
            if(vc != null) {
               variantContextsString += vc.getContig() + ":" + vc.getStart() + " ";
            } else {
               variantContextsString += "null ";
            }
        }
        System.out.println("variantContexts: " + variantContextsString);
    }
}
