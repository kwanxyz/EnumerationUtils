package com.vrtx.jdesign.gui.librarydesign.utils;

import chemaxon.calculations.clean.Cleaner;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.reaction.ConcurrentReactorProcessor;
import chemaxon.reaction.ReactionException;
import chemaxon.reaction.Reactor;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;
import chemaxon.util.iterator.MoleculeIterator;
import chemaxon.util.iterator.MoleculeIteratorFactory;
import com.vrtx.jdesign.data.librarydesign.AtomRule;
import com.vrtx.jdesign.data.librarydesign.DeprotectionGroup;
import com.vrtx.jdesign.database.JDBCElnDAO;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class EnumerationUtils {

    private static final Logger logger = LoggerFactory.getLogger(EnumerationUtils.class);


    private static List<Molecule> enumerate(RxnMolecule rxnMolecule, final Molecule[] mols, String rules, final List<DeprotectionGroup> deprotectionGroups, final boolean removeUnprotected) throws ReactionException, IOException {
        Reactor reactor = new Reactor();
        reactor.setReaction(rxnMolecule, null, "".equals(rules) ? null : rules, null, (String) null);

        ConcurrentReactorProcessor crp = new ConcurrentReactorProcessor();
        crp.setReactor(reactor);

        MoleculeIterator[] iterator = new MoleculeIterator[mols.length];
        for (int i = 0; i < mols.length; i++) {
            iterator[i] = MoleculeIteratorFactory.createMoleculeIterator(mols);
        }

        crp.setReactantIterators(iterator, ConcurrentReactorProcessor.MODE_COMBINATORIAL);

        List<Molecule> products = new ArrayList<Molecule>();
        Molecule[] reactorArray;
        while ((reactorArray = crp.react()) != null) {
            products.addAll(Arrays.asList(reactorArray));
        }

        List<Molecule> molecules = new ArrayList<>();

        for (final Molecule product : products) {
            product.dearomatize();
            Cleaner.clean(product, 2, "");

            for (int i = 0; i < product.getAtomCount(); i++) {
                MolAtom atm = product.getAtom(i);
                int idx = atm.getAtomMap();
                if (idx != 0) {
                    atm.setAtomMap(0);
                }
            }
            if(!deprotectionGroups.isEmpty()) {
                try {
                    final String origSmiles = MolExporter.exportToFormat(product, "smiles");

                    String deprotectedSmiles = JDBCElnDAO.deprotectSmiles(origSmiles, deprotectionGroups, true);
                    if(null != deprotectedSmiles){
                        molecules.add(MolImporter.importMol(deprotectedSmiles));
                    } else if(!removeUnprotected) {
                        molecules.add(product);
                    }
                } catch (Exception e1) {
                    logger.error("Error deprotecting smiles", e1);
                }
            } else {
                molecules.add(product);
            }
        }

        return molecules;
    }
}
