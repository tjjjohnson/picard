package picard.vcf.filter;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.List;

/**
 * Applies a set of hard filters to Variants and to Genotypes within a VCF.
 *
 * @author Tim Fennell
 */
public class FilterVcf extends CommandLineProgram {
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The INPUT VCF or BCF file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output VCF or BCF.")
    public File OUTPUT;

    @Option(doc="The minimum allele balance acceptable before filtering a site. Allele balance is calculated for heterozygotes as " +
            "the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote " +
            "genotypes at the same locus are measured independently. The locus is filtered if any allele balance is below the limit.")
    public double MIN_AB = 0.3d;

    @Option(doc="The minimum sequencing depth supporting a genotype before the genotype will be filtered out.")
    public int MIN_DP = 6;

    @Option(doc="The minimum genotype quality that must be achieved for a sample otherwise the genotype will be filtered out.")
    public int MIN_GQ = 20;

    @Option(doc="The maximum phred scaled fisher strand value before a site will be filtered out.")
    public double MAX_FS = 13;

    /** Constructor to default to having index creation on. */
    public FilterVcf() { this.CREATE_INDEX = true; }

    // Stock main method
    public static void main(final String[] args) {
        new FilterVcf().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final List<VariantFilter>  variantFilters = CollectionUtil.makeList(new AlleleBalanceFilter(MIN_AB), new FisherStrandFilter(MAX_FS));
        final List<GenotypeFilter> genotypeFilters = CollectionUtil.makeList(new GenotypeQualityFilter(MIN_GQ), new DepthFilter(MIN_DP));
        final VCFFileReader in = new VCFFileReader(INPUT);
        final FilterApplyingVariantIterator iterator = new FilterApplyingVariantIterator(in.iterator(), variantFilters, genotypeFilters);

        final VariantContextWriter out = new VariantContextWriterBuilder().setOutputFile(OUTPUT).build();
        final VCFHeader header = in.getFileHeader();
        header.addMetaDataLine(new VCFFilterHeaderLine("AllGtsFiltered", "Site filtered out because all genotypes are filtered out."));
        header.addMetaDataLine(new VCFFormatHeaderLine("FT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Genotype filters."));
        for (final VariantFilter filter : variantFilters) {
            for (final VCFFilterHeaderLine line : filter.headerLines()) {
                header.addMetaDataLine(line);
            }
        }

        out.writeHeader(in.getFileHeader());

        while (iterator.hasNext()) {
            out.add(iterator.next());
        }

        out.close();
        in.close();
        return 0;
    }
}
