
############################## PUBCHEM ##############################

test_that("Fetching Nicotine as name from PubChem works", {

  mol_data <- fetch_properties_from_pubchem("Nicotine")
  expect_equal(mol_data$CID, 89594)
  expect_equal(mol_data$MolecularWeight, 162.23)
})

test_that("Fetching Flufenamic acid as name from PubChem works", {

  mol_data <- fetch_properties_from_pubchem("Flufenamic acid")
  expect_equal(mol_data$CID, 3371)
  expect_equal(mol_data$MolecularWeight, 281.23)
})

test_that("Fetching Nicotine as CID from PubChem works", {

  mol_data <- fetch_properties_from_pubchem(89594)
  expect_equal(mol_data$CID, 89594)
  expect_equal(mol_data$MolecularWeight, 162.23)
})

test_that("Fetching two CIDs from PubChem works", {

  mol_data <- fetch_properties_from_pubchem(list(123, 89594))
  expect_length(mol_data, 2)
  expect_equal(mol_data$CID[1], 123)
  expect_equal(mol_data$MolecularWeight[1], 144.18)
  expect_equal(mol_data$CID[2], 89594)
  expect_equal(mol_data$MolecularWeight[2], 162.23)
})

test_that("Fetching Nicotine and  a CID from PubChem works", {

  mol_data <- fetch_properties_from_pubchem(list("Nicotine", 89594))
  expect_length(mol_data, 2)
  expect_equal(mol_data$CID[1], 89594)
  expect_equal(mol_data$MolecularWeight[1], 162.23)
  expect_equal(mol_data$CID[2], 89594)
  expect_equal(mol_data$MolecularWeight[2], 162.23)
})

test_that("Fetching Nicotine with MW and XlogP from PubChem works", {

  mol_data <- fetch_properties_from_pubchem(89594, features = c("MolecularWeight", "XLogP"))
  expect_equal(mol_data$CID, 89594)
  expect_equal(mol_data$MolecularWeight, 162.23)
  expect_equal(mol_data$XLogP, 1.2)
})

test_that("Fetching from PubChem produces a message", {

  expect_message(fetch_properties_from_pubchem(89594))
})

test_that("Fetching unknown name from PubChem fails", {

  expect_error(fetch_properties_from_pubchem("FooBardfsf"))
})

test_that("Fetching unknown CID from PubChem fails", {

  expect_error(fetch_properties_from_pubchem(-1))
})

test_that("Fetching with no name argument returns NA", {

  expect_true(is.na(fetch_properties_from_pubchem(NA)))
})

test_that("Fetching with no feature only returns CID", {

  mol_data <- fetch_properties_from_pubchem("Nicotine", NA)
  expect_length(mol_data, 1)
  expect_equal(ncol(mol_data), 1)
  expect_equal(mol_data$CID, 89594)

  mol_data <- fetch_properties_from_pubchem("Nicotine", c())
  expect_length(mol_data, 1)
  expect_equal(ncol(mol_data), 1)
  expect_equal(mol_data$CID, 89594)
})

test_that("Fetching Nicotine with unknown feature fails", {

  expect_error(fetch_properties_from_pubchem(89594, features = c("FooBar", "XLogP")))
})

############################## MOLECULE ##############################

test_that("Creating a molecule from name works", {

  mol <- molecule("Nicotine")

  expect_equal(mol$name, "Nicotine")
  expect_equal(mol$display.name, "Nicotine")
  expect_equal(mol$id, "Nicotine")
  expect_equal(mol$column.matcher,  "Nicotine")
  expect_equal(mol$pubchem.id, 89594)
  expect_equal(mol$MW, 162.23)

  expect_equal(mol$is.fraction,  F)
  expect_true(is.na(mol$fixed.unit))
})

test_that("Creating a molecule from pubchem id", {

  mol <- molecule("Test", pubchem.id = 89594)

  expect_equal(mol$name, "Test")
  expect_equal(mol$display.name, "Test")
  expect_equal(mol$id, "Test")
  expect_equal(mol$column.matcher,  "Test")
  expect_equal(mol$pubchem.id, 89594)
  expect_equal(mol$MW, 162.23)

  expect_equal(mol$is.fraction,  F)
  expect_true(is.na(mol$fixed.unit))

  expect_equal(mol$is.fraction,  F)
  expect_true(is.na(mol$fixed.unit))
})

test_that("Error handling works", {

  expect_error(molecule(""))
  expect_error(molecule(" "))
  expect_error(molecule(NA))
  expect_error(molecule(NULL))
  expect_error(molecule(12))
  expect_error(molecule(12.4))

  expect_error(molecule("Test", MW = NA))
  expect_error(molecule("Test", MW = ""))
  expect_error(molecule("Test", MW = "Test"))
  expect_error(molecule("Test", MW = c(1,2)))

  expect_error(molecule("Test", MW = 0, display.name = ""))
  expect_error(molecule("Test", MW = 0, display.name = NA))

  expect_error(molecule("Test", MW = 0, column.matcher = NA))
  expect_error(molecule("Test", MW = 0, column.matcher = ""))
  expect_error(molecule("Test", MW = 0, column.matcher = c("Test", NA)))
  expect_error(molecule("Test", MW = 0, column.matcher = c("Test", "")))

  expect_error(molecule("Test", MW = 0, id = NA))
  expect_error(molecule("Test", MW = 0, id = ""))
  expect_error(molecule("Test", MW = 0, id = c("Test", "fds")))

  expect_error(molecule("Test", MW = 0, pubchem.id = ""))
  expect_error(molecule("Test", MW = 0, pubchem.id = "Test"))
  expect_error(molecule("Test", MW = 0, pubchem.id = c(1, 2)))

  expect_error(molecule("Test", MW = 0, col = NULL))

  expect_error(molecule("Test", MW = 0, ylab = NA))
  expect_error(molecule("Test", MW = 0, ylab = NULL))

  expect_error(molecule("Test", MW = 0, lty = NULL))

  expect_error(molecule("Test", MW = 0, is.fraction = NULL))
  expect_error(molecule("Test", MW = 0, is.fraction = NA))
  expect_error(molecule("Test", MW = 0, is.fraction = 12))
  expect_error(molecule("Test", MW = 0, is.fraction = "no"))
  expect_error(molecule("Test", MW = 0, is.fraction = c(T, F)))

  expect_error(molecule("Test", MW = 0, fixed.unit = NULL))
  expect_error(molecule("Test", MW = 0, fixed.unit = 12))
  expect_error(molecule("Test", MW = 0, fixed.unit = "no"))
  expect_error(molecule("Test", MW = 0, fixed.unit = c(T, F)))
})

test_that("is_molecule works", {

  mol <- molecule("Nicotine", MW = 123)

  expect_true(is_molecule(mol))
  expect_false(is_molecule(NA))
  expect_false(is_molecule(NULL))
  expect_false(is_molecule("Hallo"))
  expect_false(is_molecule(list(mol, mol)))
})

test_that("molecule_list works", {

  mol_1 <- molecule("Nicotine", MW = 123)
  mol_2 <- molecule("Test Mol", MW = 223)
  mol_list_1 <- molecule_list(mol_1, mol_2)
  expect_true(is_molecule_list(mol_list_1))

  mol_list_2 <- molecule_list(mol_1)
  expect_true(is_molecule_list(mol_list_2))

  expect_equal(names(mol_list_1), c("Nicotine", "Test Mol"))
  expect_equal(names(mol_list_2), c("Nicotine"))

  expect_error(molecule_list())
  expect_error(molecule_list(NA))
  expect_error(molecule_list(NULL))
  expect_error(molecule_list(1))
  expect_error(molecule_list(123.4))
  expect_error(molecule_list(""))
  expect_error(molecule_list("AB"))
  expect_error(molecule_list(mol_1, ""))
  expect_error(molecule_list(mol_1, NULL))
})

test_that("molecule_ids works", {

  mol_1 <- molecule("Nicotine", MW = 123)
  mol_2 <- molecule("Test Mol", MW = 223)
  mol_list <- molecule_list(mol_1, mol_2)
  expect_equal(molecule_ids(mol_list), c("Nicotine", "Test Mol"))

  mol_list_2 <- molecule_list(mol_1)
  expect_equal(molecule_ids(mol_list_2), "Nicotine")
})

test_that("molecule_ids works", {

  mol <- molecule("Nicotine", MW = 123)

  expect_true(is_molecule_list(list(mol)))
  expect_true(is_molecule_list(list(mol, mol)))

  expect_false(is_molecule_list(list()))
  expect_false(is_molecule_list(list(mol, NA)))
  expect_false(is_molecule_list(NA))
  expect_false(is_molecule_list(NULL))
  expect_false(is_molecule_list("Hallo"))
  expect_false(is_molecule_list(list(NULL, mol)))
})

test_that("molecule_from_ids works", {

  mol_1 <- molecule("Nicotine", MW = 123)
  mol_2 <- molecule("Test Mol", MW = 223)
  mol_list <- molecule_list(mol_1, mol_2)

  one_match <- molecule_from_ids(mol_list, "Nicotine")
  expect_equal(one_match$id, "Nicotine")

  two_matches <- molecule_from_ids(mol_list, c("Nicotine", "Test Mol"))
  expect_equal(two_matches[[1]]$id, "Nicotine")
  expect_equal(two_matches[[2]]$id, "Test Mol")

  fail_two_matches <- molecule_from_ids(mol_list, c("Nicotine", "A"))
  expect_true(is.null(fail_two_matches))

  partial_two_matches <- molecule_from_ids(mol_list, c("Nicotine", "A"),
                                           partial.matches = T)
  expect_equal(partial_two_matches[[1]]$id, "Nicotine")
  expect_true(is.null(partial_two_matches[[2]]))

  partial_two_matches_2 <- molecule_from_ids(mol_list, c("Nicotine", "A"),
                                           partial.matches = T,  rm.na = T)
  expect_equal(partial_two_matches_2$id, "Nicotine")
})

test_that("has_molecules works", {

  mol_1 <- molecule("Nicotine", MW = 123)
  mol_2 <- molecule("Test Mol", MW = 223)
  mol_list <- molecule_list(mol_1, mol_2)

  expect_true(has_molecules(mol_list, c("Nicotine", "Test Mol")))
  expect_true(has_molecules(mol_list, c("Nicotine")))
  expect_true(has_molecules(mol_list, c("Test Mol")))

  expect_false(has_molecules(mol_list, c("Nicotine", "A")))
  expect_false(has_molecules(mol_list, c("A")))
  expect_false(has_molecules(mol_list, c("")))
})

test_that("merge_molecule_lists works", {

  mol_1 <- molecule("Nicotine", MW = 123)
  mol_2 <- molecule("Test Mol", MW = 223)
  mol_list_1 <- molecule_list(mol_1, mol_2)

  mol_3 <- molecule("Nicotine", MW = 123)
  mol_4 <- molecule("Test Mol 2", MW = 223)
  mol_list_2 <- molecule_list(mol_3, mol_4)

  self_merge <- merge_molecule_lists(mol_list_1, mol_list_1)
  expect_length(self_merge, 2)
  expect_true(has_molecules(self_merge, c("Nicotine", "Test Mol")))

  merged <- merge_molecule_lists(mol_list_1, mol_list_2)
  expect_length(merged, 3)
  expect_true(has_molecules(merged, c("Nicotine", "Test Mol", "Test Mol 2")))
})

test_that("is_unique_molecule_list works", {

  mol_1 <- molecule("Nicotine", MW = 123)
  mol_2 <- molecule("Test Mol", MW = 223)
  mol_list_1 <- molecule_list(mol_1, mol_2)

  expect_true(is_unique_molecule_list(mol_list_1))
  expect_true(is_unique_molecule_list(molecule_list(mol_1)))

  mol_list_2 <- molecule_list(mol_1, mol_2, mol_1)
  expect_false(is_unique_molecule_list(mol_list_2))
})


############################## MASTER ##############################

test_that("creating a master object works", {

  mols <- molecule_list(molecule("Test", MW = 100))

  m1 <- master("first", molecules = mols)
  expect_equal(m1$id, "first")
  expect_true(is.na(m1$reference))
  expect_length(m1$group.names, 0)
  expect_equal(m1$molecules[[1]]$id, "Test")

  m2 <- master("second", "ref",  groups = 2, molecules = mols)
  expect_equal(m2$id, "second")
  expect_equal(m2$reference, "ref")
  expect_length(m2$group.names, 2)

  m3 <- master("second", "ref",  groups = 2, molecules = molecule("Test", MW = 100))
  expect_equal(m3$molecules[[1]]$id, "Test")
})

test_that("creating a master object is robust", {

  expect_error(master("No molecule", NA))

  dummy_mol <- molecule("Test", MW = 100)
  expect_error(master("", molecule_list(dummy_mol, dummy_mol)))

  expect_error(master("", dummy_mol))
  expect_error(master(" ", dummy_mol))
  expect_error(master(NA, dummy_mol))

  expect_error(master("test", groups = -1, dummy_mol))
  expect_error(master("test", groups = "test", dummy_mol))

  expect_error(master("test", groups = 2, molecules = NULL))
  expect_error(master("test", groups = 2, molecules = NA))
  expect_error(master("test", groups = 2, molecules = "Hallo"))
})

test_that("add_master_entry works - basic example", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m)
  expect_equal(nrow(m$data), 1)

  m <- add_master_entry(m)
  expect_equal(nrow(m$data), 2)

  expect_true(all(is.na(m$data[1,])))
  expect_true(all(is.na(m$data[2,])))
})

test_that("add_master_entry works - basic error", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))

  expect_error(add_master_entry("Hallo"))
  expect_error(add_master_entry(NA))
  expect_error(add_master_entry(NULL))
})

test_that("add_master_entry works - add pop", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m, pop.file = "Pop", pop.molecules = mol_1)
  m <- add_master_entry(m, pop.file = "Pop 2", pop.molecules = mol_2)
  m <- add_master_entry(m, pop.file = "Pop 3", pop.molecules = molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 3)
})

test_that("add_master_entry works - add pop error", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  expect_error(add_master_entry(m, pop.id = "Pop", pop.molecules = molecule_list(mol_1, mol_1)))
  expect_error(add_master_entry(m, pop.id = "Pop", pop.molecules = NA))
  expect_error(add_master_entry(m, pop.pop.idname = "Pop", pop.molecules = NULL))
  expect_error(add_master_entry(m, pop.id = "Pop", pop.molecules = molecule("Too", MW = 100)))
  expect_error(add_master_entry(m, pop.id = NA, pop.molecules = molecule_list(mol_1, mol_2)))
  expect_error(add_master_entry(m, pop.id = NULL, pop.molecules = molecule_list(mol_1, mol_2)))
  expect_error(add_master_entry(m, pop.id = c("1", "2"), pop.molecules = molecule_list(mol_1, mol_2)))
  expect_error(add_master_entry(m, pop.id = list("1", "2"), pop.molecules = molecule_list(mol_1, mol_2)))
})

test_that("add_master_entry works - add sim", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m, sim.file = "Sim", sim.molecules = mol_1)
  m <- add_master_entry(m, sim.file = "Sim 2", sim.molecules = mol_2)
  m <- add_master_entry(m, sim.file = "Sim 3", sim.molecules = molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 3)
})

test_that("add_master_entry works - add sim error", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  expect_error(add_master_entry(m, sim.file = "Sim", sim.molecules = molecule_list(mol_1, mol_1)))
  expect_error(add_master_entry(m, sim.file = "Sim", sim.molecules = NA))
  expect_error(add_master_entry(m, sim.file = "Sim", sim.molecules = NULL))
  expect_error(add_master_entry(m, sim.file = "Sim", sim.molecules = molecule("Too", MW = 100)))
  expect_error(add_master_entry(m, sim.file = NA, sim.molecules = molecule_list(mol_1, mol_2)))
  expect_error(add_master_entry(m, sim.file = NULL, sim.molecules = molecule_list(mol_1, mol_2)))
  expect_error(add_master_entry(m, sim.file = c("1", "2"), sim.molecules = molecule_list(mol_1, mol_2)))
  expect_error(add_master_entry(m, sim.file = list("1", "2"), sim.molecules = molecule_list(mol_1, mol_2)))
})

test_that("add_master_entry works - add obs", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m, obs.ids = 1)
  m <- add_master_entry(m, obs.ids = "A")
  m <- add_master_entry(m, obs.ids = list(1, "A"))
  m <- add_master_entry(m, obs.ids = c("B", "A"))
  expect_equal(nrow(m$data), 4)
})

test_that("add_master_entry works - add obs error", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  expect_error(add_master_entry(m, obs.ids = NULL))
  expect_error(add_master_entry(m, obs.ids = c(1, NA)))
  expect_error(add_master_entry(m, obs.ids = c(NA, "Test")))
  expect_error(add_master_entry(m, obs.ids = list(NA, "Test")))
})

test_that("add_master_entry works - add header", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m, plot.header = "")
  m <- add_master_entry(m, plot.header = "Test")
  m <- add_master_entry(m, plot.header = 12)
  expect_equal(nrow(m$data), 3)
  expect_equal(m$data$Plot_Header, c("", "Test", "12"))
})

test_that("add_master_entry works - add header error", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  expect_error(add_master_entry(m, plot.header = NULL))
  expect_error(add_master_entry(m, plot.header = c(1,2)))
  expect_error(add_master_entry(m, plot.header = c("A","B")))
})

test_that("add_master_entry works - x_range", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m, plot.x.range = NA)
  m <- add_master_entry(m, plot.x.range = units::as_units(c(1,2), "h"))
  m <- add_master_entry(m, plot.x.range = units::as_units(c(NA, 2), "h"))
  m <- add_master_entry(m, plot.x.range = units::as_units(c(1, NA), "h"))
  expect_equal(nrow(m$data), 4)
  expect_equal(m$data$Plot_X_Min, c(NA, 1, NA, 1))
  expect_equal(m$data$Plot_X_Max, c(NA, 2, 2, NA))
  expect_equal(m$data$Plot_X_Unit, c(NA, "h", "h", "h"))
})

test_that("add_master_entry works - x_range errors", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  expect_error(add_master_entry(m, plot.x.range = c(1,2)))
  expect_error(add_master_entry(m, plot.x.range = units::as_units(c(1,2,3), "h")))
})

test_that("add_master_entry works - y_range", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m, plot.y.range = NA)
  m <- add_master_entry(m, plot.y.range = units::as_units(c(1,2), "h"))
  m <- add_master_entry(m, plot.y.range = units::as_units(c(NA, 2), "h"))
  m <- add_master_entry(m, plot.y.range = units::as_units(c(1, NA), "h"))
  expect_equal(nrow(m$data), 4)
  expect_equal(m$data$Plot_Y_Min, c(NA, 1, NA, 1))
  expect_equal(m$data$Plot_Y_Max, c(NA, 2, 2, NA))
  expect_equal(m$data$Plot_Y_Unit, c(NA, "h", "h", "h"))
})

test_that("add_master_entry works - y_range errors", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2))
  expect_equal(nrow(m$data), 0)

  expect_error(add_master_entry(m, plot.y.range = c(1,2)))
  expect_error(add_master_entry(m, plot.y.range = units::as_units(c(1,2,3), "h")))
})

test_that("add_master_entry works - groups", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2), groups = 2)
  expect_equal(nrow(m$data), 0)

  m <- add_master_entry(m, groups = NA)
  m <- add_master_entry(m, groups = list(NA, "Bar"))
  m <- add_master_entry(m, groups = c("Foo", NA))
  m <- add_master_entry(m, groups = c("Hey", 1))

  expect_equal(nrow(m$data), 4)
  expect_equal(m$data$Group_1, c(NA, NA, "Foo", "Hey"))
  expect_equal(m$data$Group_2, c(NA, "Bar", NA, "1"))
})

test_that("add_master_entry works - groups errors", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)

  m <- master("Test", molecule_list(mol_1, mol_2), groups = 2)
  expect_equal(nrow(m$data), 0)

  expect_error(add_master_entry(m, groups = c(1, 2, 3)))
})

test_that("bind_master works", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)
  mol_3 <- molecule("Baz", MW = 300)

  m_1 <- master("one", molecules = molecule_list(mol_1, mol_2), groups = 3)
  m_1 <- add_master_entry(m_1, pop.file = "Pop", pop.molecules = mol_1, groups = c("A", NA, "B"))

  m_2 <- master("one", molecules = molecule_list(mol_1, mol_3), groups = 2)
  m_2 <- add_master_entry(m_2, sim.file = "Sim", sim.molecules = mol_3, groups = c(NA, "C"))

  merged <- bind_master(m_1, m_2)

  expect_equal(nrow(merged$data), 2)
  expect_length(merged$molecules, 3)
  expect_equal(molecule_ids(merged$molecules), c("Foo", "Bar", "Baz"))
  expect_length(merged$group.names, 3)
  expect_equal(merged$data$Pop_File, c("Pop", NA))
  expect_equal(merged$data$Pop_Id, c("Pop", NA))
  expect_equal(merged$data$Sim_File, c(NA, "Sim"))
  expect_equal(merged$data$Sim_Id, c(NA, "Sim"))
  expect_equal(merged$data$Pop_Mols, c("Foo", NA))
  expect_equal(merged$data$Sim_Mols, c(NA, "Baz"))
  expect_equal(merged$data$Group_1, c("A", NA))
  expect_equal(merged$data$Group_2, c(NA, "C"))
  expect_equal(merged$data$Group_3, c("B", NA))
})


test_that("bind_master works with errors", {

  mol_1 <- molecule("Foo", MW = 100)
  mol_2 <- molecule("Bar", MW = 200)
  mol_3 <- molecule("Baz", MW = 300)

  m_1 <- master("one", molecules = molecule_list(mol_1, mol_2), groups = 3)
  m_1 <- add_master_entry(m_1, pop.file = "Pop", pop.molecules = mol_1, groups = c("A", NA, "B"))

  m_2 <- master("one", molecules = molecule_list(mol_1, mol_3), groups = 2)
  m_2 <- add_master_entry(m_2, sim.file = "Sim", sim.molecules = mol_3, groups = c(NA, "C"))

  expect_error(bind_master(m_2, m_3))
})


test_that("groups works", {

  mol <- molecule("Mol 1", MW = 100)
  m1 <- master("test", mol)
  m2 <- master("test 2", mol, "my ref", groups = 3)
  expect_true(is.null(groups(m1)))
  expect_equal(groups(m2), c("Group_1", "Group_2", "Group_3"))
})


test_that("entry works - fn call", {

  mol <- molecule("Mol 1", MW = 100)
  m1 <- master("test", mol)
  m1 <- add_master_entry(m1, "Pop 1", pop.molecules = mol)
  m1 <- add_master_entry(m1, "Pop 2", pop.molecules = mol)
  m1 <- add_master_entry(m1, "Pop 3", pop.molecules = mol)
  m1 <- add_master_entry(m1, "Pop 4", pop.molecules = mol)


  m_filtered <- entry(m1, 1)
  expect_equal(m_filtered$Pop_Id, "Pop 1")

  m_filtered_2 <- entry(m1, c(1,2))
  expect_equal(m_filtered_2[[1]]$Pop_Id, "Pop 1")
  expect_equal(m_filtered_2[[2]]$Pop_Id, "Pop 2")
})


test_that("entry works - operator call", {

  mol <- molecule("Mol 1", MW = 100)
  m1 <- master("test", mol)
  m1 <- add_master_entry(m1, "Pop 1", pop.molecules = mol)
  m1 <- add_master_entry(m1, "Pop 2", pop.molecules = mol)
  m1 <- add_master_entry(m1, "Pop 3", pop.molecules = mol)
  m1 <- add_master_entry(m1, "Pop 4", pop.molecules = mol)


  m_filtered <- m1[1]
  expect_equal(m_filtered$Pop_Id, "Pop 1")

  m_filtered_2 <- m1[1,2]
  expect_equal(m_filtered_2[[1]]$Pop_Id, "Pop 1")
  expect_equal(m_filtered_2[[2]]$Pop_Id, "Pop 2")
})

