package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path"
	"path/filepath"
	"strconv"
	"strings"
)

type Cgene struct {
	StartBP       int64
	EndBP         int64
	Strand        string
	ProteinLength int64
	PID           string
	Name          string
	Product       string
}

type Organism struct {

	// ------------- organism info -------------------------------------------------
	OrganismName string

	// ------------- protein(gene) info --------------------------------------------
	PttFilename string
	Ngenes      int64
	Ptt         []Cgene

	// ------------- DNA info ---------------------------------------
	DNAFilename  string
	Nnucleotides int64
	HeaderLine   string
	DNA          string
}

func getDataPath() string {
	ex, err := os.Executable()
	if err != nil {
		panic(err)
	}
	excuteDir := filepath.Dir(ex)
	return path.Join(excuteDir, "data")
}

func (organism *Organism) LoadProteinTable(organismName string) bool {
	log.Printf("\tAttempt to Load Protein Table (list of genes) for organism:%s\n", organismName)
	organism.PttFilename = path.Join(getDataPath(), fmt.Sprintf("%s.ptt", organismName))
	log.Printf("\tAttempt to open filename: %s\n", organism.PttFilename)
	if strings.TrimSpace(organism.PttFilename) == "" {
		log.Fatalln("[LoadProteinTable]:PttFilename can't be null")
		return false
	}
	if _, err := os.Stat(organism.PttFilename); os.IsNotExist(err) {
		log.Fatalln("[LoadProteinTable]:PttFilename not exists")
		return false
	}
	pttBytes, err := ioutil.ReadFile(organism.PttFilename)
	if err != nil {
		log.Fatalln("[LoadProteinTable]:", err)
		return false
	}

	ptt_text := string(pttBytes)

	line_texts := strings.Split(strings.TrimSpace(ptt_text), "\n")

	organism.Ptt = make([]Cgene, 0)

	for idx, line := range line_texts {
		if idx == 0 {
			// parse header
			continue
		}
		if idx == 1 {
			// proteins length
			protein_len_texts := strings.Split(line, " ")
			if len(protein_len_texts) == 0 {
				log.Println("[LoadProteinTable]：parse protein length error")
				continue
			}

			proteion_len, err := strconv.ParseInt(protein_len_texts[0], 10, 64)
			if err != nil {
				log.Println("[LoadProteinTable]：parse protein length error ", err)
				continue
			}
			organism.Ngenes = proteion_len
			log.Printf("\tNumber of lines (excl. header) in %s.ptt is: %d\n", organismName, organism.Ngenes)
			continue
		}
		if idx == 2 {
			// column header
			continue
		}

		// log.Println(line)
		gene_line_texts := strings.Split(line, "\t")
		if len(gene_line_texts) < 9 {
			log.Println("[LoadProteinTable]: data error")
			continue
		}

		protein_len, err := strconv.ParseInt(gene_line_texts[2], 10, 64)

		if err != nil {
			log.Println("[LoadProteinTable]: data error", err)
			continue
		}

		location_texts := strings.Split(gene_line_texts[0], "..")
		start_bp, err := strconv.ParseInt(location_texts[0], 10, 64)
		end_bp, err := strconv.ParseInt(location_texts[1], 10, 64)

		gene := Cgene{
			StartBP:       start_bp,
			EndBP:         end_bp,
			Strand:        gene_line_texts[1],
			ProteinLength: protein_len,
			PID:           gene_line_texts[3],
			Name:          gene_line_texts[4],
			Product:       gene_line_texts[8],
		}

		organism.Ptt = append(organism.Ptt, gene)
	}
	log.Printf("\tProtein table %s.ptt successfully opened and loaded.\n\n", organismName)
	return true
}

func (organism *Organism) printProteinTable2stdout() {

}

func (organism *Organism) PrintProteinTable2stdout(n int64) {

}

func (organism *Organism) LoadDNA(organismName string) bool {
	log.Printf("\tAttempt to Load DNA for organism:%s\n", organismName)
	organism.DNAFilename = path.Join(getDataPath(), fmt.Sprintf("%s.fna", organismName))
	log.Printf("\tAttempt to open filename: %s\n", organism.DNAFilename)
	if strings.TrimSpace(organism.DNAFilename) == "" {
		log.Fatalln("[LoadDNA]:DNAFilename can't be null")
		return false
	}
	if _, err := os.Stat(organism.DNAFilename); os.IsNotExist(err) {
		log.Fatalln("[LoadDNA]:DNAFilename not exists")
		return false
	}
	dnaBytes, err := ioutil.ReadFile(organism.DNAFilename)
	if err != nil {
		log.Fatalln("[LoadDNA]:", err)
		return false
	}

	dna_text := string(dnaBytes)
	dna_line_texts := strings.SplitN(dna_text, "\n", 2)

	if len(dna_line_texts) < 2 {
		log.Fatalln("[LoadDNA]: parse failed")
		return false
	}

	organism.HeaderLine = dna_line_texts[0]

	header_infos := strings.Split(organism.HeaderLine, "|")

	if len(header_infos) > 1 {
		organism.Nnucleotides, err = strconv.ParseInt(header_infos[1], 10, 64)
		log.Printf("\tNumber of lines of DNA (excl. header) in %s.fna is: %d\n", organismName, organism.Nnucleotides)
	}

	organism.DNA = strings.ReplaceAll(dna_line_texts[1], "\n", "")

	log.Printf("\tDNA %s.fna successfully opened and loaded.\n\n", organismName)
	return true
}

func (organism *Organism) headDNA2stdout(n int64) {

}

func (organism *Organism) bsearch(search string) {

	dna_search_index := int64(strings.Index(organism.DNA, search))

	dna_search_end_index := dna_search_index + int64(len(search))

	ptts := make([]Cgene, 0)

	for ptt_idx, ptt := range organism.Ptt {
		if ptt.StartBP <= dna_search_index && ptt.EndBP >= dna_search_end_index {
			ptts = append(ptts, ptt)
			break
		}
		if ptt.StartBP > dna_search_end_index {
			if ptt_idx > 0 {
				last_ptt := organism.Ptt[ptt_idx-1]
				if last_ptt.StartBP <= dna_search_index {
					ptts = append(ptts, last_ptt)
					ptts = append(ptts, ptt)
					break
				}
			}

		}
	}

	if len(ptts) == 0 {
		log.Printf("Not Found\n")
	} else if len(ptts) == 1 {
		ptt := ptts[0]
		log.Printf("%s found at bp %d\n", search, dna_search_index)
		log.Printf("WITHIN GENE Name: %s Start: %d End: %d PID: %s\n", ptt.Name, ptt.StartBP, ptt.EndBP, ptt.PID)
	} else if len(ptts) == 2 {
		up := ptts[0]
		down := ptts[1]
		log.Printf("\t%s found at bp %d\n", search, dna_search_index)
		log.Printf("\t\tIn between genes:\n")
		log.Printf("\t\t\tUPSTREAM Name: %s Start: %d End: %d PID: %s\n", up.Name, up.StartBP, up.EndBP, up.PID)
		log.Printf("\t\t\tand\n")
		log.Printf("\t\t\tDOWNSTREAM Name: %s Start: %d End: %d PID: %s\n", down.Name, down.StartBP, down.EndBP, down.PID)
	}

}

func printLine() {
	log.Println("----------------------------------------------------------")
}

func printDoubleLine() {
	log.Println("==========================================================================")
}

func main() {
	log.SetFlags(0)
	log.SetOutput(new(logWriter))

	for {
		printDoubleLine()

		fmt.Printf("Enter the Organism name (e.g., Ecoli): ")

		var orgaismName string
		fmt.Scanln(&orgaismName)

		if len(orgaismName) == 0 {
			continue
		}

		organism := &Organism{}

		organism.LoadProteinTable(orgaismName)
		organism.LoadDNA(orgaismName)

		for {
			fmt.Printf("What motif would you like to find >> ")
			var motif string
			fmt.Scanln(&motif)
			organism.bsearch(motif)
			printLine()

			var response string
			var breakParent bool = false

			for len(response) == 0 {
				fmt.Printf("Do you want to search for another motif? (Y/n): ")

				fmt.Scanln(&response)

				printLine()

				if len(response) > 0 {
					if response[0] == 'N' || response[0] == 'n' {
						breakParent = true
						break
					}
				}
			}
			if breakParent {
				break
			}
		}

		fmt.Printf("Do you want to try another organism? (Y/n):")
		var response string
		fmt.Scanln(&response)

		printLine()
		if len(response) > 0 {
			if response[0] == 'N' || response[0] == 'n' {
				break
			}
		}

	}

	printDoubleLine()
}

type logWriter struct {
}

func (writer logWriter) Write(bytes []byte) (int, error) {
	return fmt.Print(string(bytes))
}
