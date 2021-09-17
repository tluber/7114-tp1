import java.io.BufferedWriter
import java.io.File

fun main(args: Array<String>) {

    parseFile()
    process()
    printResults()
}

private val sortedGarments = mutableListOf<Garment>()
private val results = mutableListOf<Pair<Int, Int>>()

private fun parseFile() {
    val fileRoute = "src/main/resources/problem"
    val lines: List<String> = File(fileRoute).readLines()

    val garmentIncompatibilities = lines.filter { line -> line.startsWith('e') }.map { line -> line.split(" ") }
    val washTimes = lines.filter { line -> line.startsWith('n') }.map { line -> line.split(" ") }

    val garments = mutableListOf<Garment>()

    for (washTime in washTimes) {
        val id = washTime[1].toInt()
        val time = washTime[2].toInt()
        val incompatibilities =
            garmentIncompatibilities.filter { list -> list[1] == washTime[1] }.map { item -> item[2].toInt() }

        val garment = Garment(id, time, incompatibilities)
        garments.add(garment)
    }
    sortedGarments.addAll(garments.sortedBy { it.washTime })
}

private fun process() {
    var n = 1
    for (garment in sortedGarments) {
        if (!isInResults(garment.id)) {
            results.addAll(getPairs(garment.id, n))
            n++
        }
    }
}

private fun getPairs(garmentId: Int, washId: Int): List<Pair<Int, Int>> {
    val result = mutableListOf<Pair<Int, Int>>()
    val list = mutableListOf<Garment>()

    list.addAll(sortedGarments.filter { it.id == garmentId })
    for (garment in sortedGarments) {
        if (!isInResults(garment.id)) {
            if (!list.map { it.isCompatible(garment.id) }.any { !it }) {
                result.add(Pair(garment.id, washId))
                list.add(garment)
            }
        }
    }

    return result
}

private fun isInResults(garmentId: Int): Boolean {
    for (result in results) {
        if (result.first == garmentId) {
            return true
        }
    }
    return false
}

private fun printResults() {
    val fileName = "src/main/resources/results.txt"

    val file = File(fileName)

    val isFileCreated: Boolean = file.createNewFile()

    if (isFileCreated) {
        file.bufferedWriter().use { out ->
            results.forEach {
                out.writeLn("${it.first} ${it.second}")
            }
        }
    } else {
        println("Results.txt already exists.")
    }
}

data class Comment(val line: String)
data class Problem(val comment: Comment, val garments: Int, val incompatibilities: Int)
data class Garment(val id: Int, val washTime: Int, val incompatibilities: List<Int>) {

    fun isCompatible(id: Int): Boolean {
        return !this.incompatibilities.contains(id)
    }
}

fun BufferedWriter.writeLn(line: String) {
    this.write(line)
    this.newLine()
}