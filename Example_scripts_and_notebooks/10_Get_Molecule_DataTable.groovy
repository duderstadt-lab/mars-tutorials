// This script retrieves the molecule DataTable for a specified UUID

#@ MoleculeArchive archive
#@ String UID
#@ output MarsTable datatable

datatable = archive.get(UID).getDataTable().clone()
