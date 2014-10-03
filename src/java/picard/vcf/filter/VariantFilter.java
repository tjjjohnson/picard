package picard.vcf.filter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;

import java.util.List;

/**
 * Interface for classes that can generate filters for VariantContexts. The contract is that a
 * VariantContext is provided, and if the variant should be filtered out then the filter string
 * should be returned, otherwise null.
 *
 * @author tfennell
 */
public interface VariantFilter {
    public String filter(final VariantContext ctx);

    public List<VCFFilterHeaderLine> headerLines();
}
