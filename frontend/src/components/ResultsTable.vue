<template>
  <table border="1">
    <thead>
      <tr>
        <th v-for="(value, key) in dictionary" :key="key">{{ key }}</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td v-for="(value, key) in dictionary" :key="key">{{ value }}</td>
      </tr>
    </tbody>
  </table>
</template>

<script>
export default {
  props: {
    dictionary: {
      type: Object,
      required: true
    }
  }
}
</script>
