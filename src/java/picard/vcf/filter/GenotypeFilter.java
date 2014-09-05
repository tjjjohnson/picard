package picard.vcf.filter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Map;

/**
 *
 */
public interface GenotypeFilter {
    public String filter(final VariantContext ctx, final Genotype gt);
}
