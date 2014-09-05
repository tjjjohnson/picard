package picard.vcf.filter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import static htsjdk.samtools.util.CollectionUtil.MultiMap;

import htsjdk.samtools.util.ListMap;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * Iterator that dynamically applies filter strings to VariantContext records supplied by an underlying
 * iterator.  Returns all records from the underlying stream and does not remove any.
 *
 * @author tfennell
 */
public class FilterApplyingVariantIterator implements CloseableIterator<VariantContext> {
    private final Iterator<VariantContext> iterator;
    private final VariantFilter[] filters;
    private final GenotypeFilter[] gtFilters;


    public FilterApplyingVariantIterator(final Iterator<VariantContext> iterator,
                                         final Collection<VariantFilter> filters,
                                         final Collection<GenotypeFilter> gtFilters) {
        this.iterator = iterator;
        this.filters = filters.toArray(new VariantFilter[filters.size()]);
        this.gtFilters = gtFilters.toArray(new GenotypeFilter[gtFilters.size()]);
    }

    /**
     * Provides the next record from the underlying iterator after applying filter strings generated
     * by the set of filters in use by the iterator.
     */
    @Override
    public VariantContext next() {
        final VariantContext ctx = this.iterator.next();
        final Set<String> filterStrings = new HashSet<String>();

        // Collect variant level filters
        for (final VariantFilter filter : this.filters) {
            final String val = filter.filter(ctx);
            if (val != null) filterStrings.add(val);
        }

        // Collect genotype level filters in a Map of Sample -> List<filter string>
        final ListMap<String,String> gtFilterStrings = new ListMap<String,String>();
        for (final Genotype gt : ctx.getGenotypes()) {
            for (final GenotypeFilter filter : gtFilters) {
                gtFilterStrings.add(gt.getSampleName(), filter.filter(ctx, gt));
            }
        }

        // If we haven't accumulated any filters at all, just pass the record through
        if (filterStrings.isEmpty() && gtFilterStrings.isEmpty()) return ctx;

        final VariantContextBuilder builder = new VariantContextBuilder(ctx);
        if (!filterStrings.isEmpty()) builder.filters(filterStrings);

        if (!gtFilterStrings.isEmpty()) {
            // Apply filters to the necessary genotypes
            builder.noGenotypes();
            final List<Genotype> newGenotypes = new ArrayList<Genotype>(ctx.getNSamples());
            for (final Genotype gt : ctx.getGenotypes()) {
                final List<String> filters = gtFilterStrings.get(gt.getSampleName());
                if (filters == null || filters.isEmpty()) {
                    newGenotypes.add(gt);
                }
                else {
                    final GenotypeBuilder gtBuilder = new GenotypeBuilder(gt);
                    gtBuilder.filters(filters);
                    newGenotypes.add(gtBuilder.make());
                }
            }
            builder.genotypes(newGenotypes);


            // If all genotypes are filtered apply a site level filter
            if (gtFilterStrings.size() == ctx.getNSamples()) {
                builder.filter("AllGtsFiltered");
            }
        }

        return builder.make();
    }

    @Override public boolean hasNext() { return this.iterator.hasNext(); }
    @Override public void close() { CloserUtil.close(this.iterator); }
    @Override public void remove() { throw new UnsupportedOperationException("remove() not supported by FilterApplyingVariantIterator."); }
}
