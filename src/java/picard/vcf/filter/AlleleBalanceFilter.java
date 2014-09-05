package picard.vcf.filter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Filters out a record if the allele balance is out of a defined range across all samples.  In the case of Heterozygous
 * genotypes the threshold is set as the minimum fraction of the data drawn from the less-represented allele - e.g. 0.3 would
 * set that whichever allele has lower representation across all heterozygous individuals must account for at least 30% of the
 * total observations.  For homozygous genotypes the number is the *upper* bound on the lower represented allele.
 */
public class AlleleBalanceFilter implements VariantFilter {
    private final int hetAlleleBalance, homAlleleBalance;

    public AlleleBalanceFilter(final int hetAlleleBalance, final int homAlleleBalance) {
        this.hetAlleleBalance = hetAlleleBalance;
        this.homAlleleBalance = homAlleleBalance;
    }

    private static class Counts { int allele1; int allele2; }

    @Override
    public String filter(final VariantContext ctx) {
        final Map<List<Allele>, Counts> countsMap = new HashMap<List<Allele>, Counts>();

        final Map<Allele, Integer> alleleIndices = new HashMap<Allele,Integer>();
//        for (Allele allele : ctx.getAlleles()) all

        for (final Genotype gt : ctx.getGenotypes()) {
            final List<Allele> alleles = gt.getAlleles();
            if (gt.isHet()) {
                Counts counts = countsMap.get(alleles);
                if (counts == null) {
                    counts = new Counts();
                    countsMap.put(alleles, counts);
                }




            }
        }



        return null;
    }
}

