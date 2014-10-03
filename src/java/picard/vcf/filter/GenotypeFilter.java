package picard.vcf.filter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Map;

/**
 * An interface for classes that perform Genotype filtration. Implementations are expected to take in a VariatnContext
 * and a single Genotype and return either null (for no filter) or a specific filter string.
 *
 * @author Tim Fennell
 */
public interface GenotypeFilter {
    public String filter(final VariantContext ctx, final Genotype gt);
}
