<template>
  <div>
    <h2>Home Page</h2>
    <input type='file' @change='onFileChange' />
    <button @click='upload'>Upload File</button>
    <button v-if="resFileLink" @click="() => downloadResults(resFileLink, resFile)">Download Results</button>
    <div v-if='waiting'>Processing...</div>
    <div v-if="result !== null" class="centered-container">
        <ResultsTable :dictionary="result" />
    </div>
    <div v-if="result !== null">
        <label for="steps">Steps:</label>
        <input type="number" id="steps" v-model.number="steps" />
        <label for="fmax">Fmax:</label>
        <input type="number" id="fmax" step="0.01" v-model.number="fmax" />
        <button @click="relax">Perform Relax</button>
        <button v-if="trajFileLink" @click="() => downloadResults(trajFileLink, trajFile)">Download Trajectory</button>
        <button v-if="relaxFileLink" @click="() => downloadResults(relaxFileLink, relaxFile)">Download Relaxed Structure</button>
    </div>
    <div v-if='waiting_relax'>Relaxing...</div>
    <div v-if="result_relax !== null" class="centered-container">
        <ResultsTable :dictionary="result_relax" />
    </div>
  </div>
</template>

<style>
.centered-container {
  display: flex;
  justify-content: center;
  align-items: center;
  min-height: 100%; /* Optional: Set a minimum height if needed */
}
</style>

<script>

import axios from 'axios'
import ResultsTable from '@/components/ResultsTable.vue'

export default {
  components: {
    ResultsTable
  },
  data () {
    return {
      file: null,
      resFile: null,
      resFileLink: '',
      waiting: false,
      result: null,
      steps: 100,
      fmax: 0.1,
      trajFile: null,
      trajFileLink: '',
      relaxFile: null,
      relaxFileLink: '',
      waiting_relax: false,
      result_relax: null
    }
  },
  methods: {
    onFileChange (e) {
      this.file = e.target.files[0]
      this.resFile = null
      this.resFileLink = ''
      this.waiting = false
      this.result = null
      this.trajFile = null
      this.trajFileLink = ''
      this.relaxFile = null
      this.relaxFileLink = ''
      this.waiting_relax = false
      this.result_relax = null
    },
    async upload () {
      if (!this.file) {
        alert('Please select a file for upload')
        return
      }

      const formData = new FormData()
      formData.append('file', this.file)

      try {
        this.waiting = true
        const response = await fetch('http://localhost:5000/pes_calc', {
          method: 'POST',
          body: formData
        })
        const data = await response.json()
        if (data.status === 'success') {
          this.resFileLink = `http://localhost:5000/download/${data.res_file}`
          this.result = data.result
          this.resFile = data.res_file
        } else {
          alert(data.message)
        }
      } catch (error) {
        alert('Error in uploading file')
      } finally {
        this.waiting = false
      }
    },
    async relax () {
      if (!this.file) {
        alert('Please select a file for upload')
        return
      }

      const formData = new FormData()
      formData.append('file', this.file)
      formData.append('steps', this.steps)
      formData.append('fmax', this.fmax)

      try {
        this.waiting_relax = true
        const response = await fetch('http://localhost:5000/relax', {
          method: 'POST',
          body: formData
        })
        const data = await response.json()
        if (data.status === 'success') {
          this.relaxFileLink = `http://localhost:5000/download/${data.relaxed_file}`
          this.relaxFile = data.relaxed_file
          this.trajFileLink = `http://localhost:5000/download/${data.traj_file}`
          this.trajFile = data.traj_file
          this.result_relax = data.result
        } else {
          alert(data.message)
        }
      } catch (error) {
        alert('Error in uploading file')
      } finally {
        this.waiting_relax = false
      }
    },
    async downloadResults (downloadLink, fileName) {
      try {
        const response = await axios.get(downloadLink, { responseType: 'blob' })
        const url = window.URL.createObjectURL(new Blob([response.data]))
        const link = document.createElement('a')
        link.href = url
        link.setAttribute('download', fileName) // Use the original file name for the downloaded file
        document.body.appendChild(link)
        link.click()
        link.remove()
      } catch (error) {
        alert('Error downloading the file')
      }
    }
  }
}
</script>
